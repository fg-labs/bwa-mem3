/* bwa_shm_stage holds the stager lock only for registry bookkeeping, marks a
 * segment complete only after the pack, and recovers from an interrupted or
 * failed stage. Every check stages a small index under private per-run
 * registry, lock and segment names, so the host registry is never touched.
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
#include "framework/temp_dir.h"

#include <algorithm>
#include <semaphore.h>
#include <signal.h>
#include <sys/ipc.h>
#include <sys/mman.h>
#include <sys/shm.h>
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
// g_fake_frsize is the allocation unit the probe reports (f_bavail counts
// whole units of it).
static int g_wait_s = 30;
static long long g_fake_free_bytes = -1;
static unsigned long g_fake_frsize = 1;
static int fake_statvfs(const char *, struct statvfs *v) {
    if (g_fake_free_bytes < 0) return -1;
    memset(v, 0, sizeof *v);
    v->f_frsize = g_fake_frsize;
    v->f_bavail = (fsblkcnt_t)((unsigned long long)g_fake_free_bytes / g_fake_frsize);
    return 0;
}
// How long `shm -d` waits for a stager lock that another process holds.
static int g_destroy_wait_ms = 5000;
#define BWA_SHM_DESTROY_LOCK_WAIT_MS g_destroy_wait_ms
// A 1 GB stage is test-controlled too: by default it makes an ordinary System V
// segment of whole 64 KiB "pages", which needs no huge page pool, and
// g_token != 0 fixes the stage token, and so the segment's key.
static int g_huge_flags = 0;
static unsigned long long g_huge_page = 65536;
static uint64_t g_token = 0;
static uint64_t stage_token(void);
static uint64_t test_token(void) { return g_token != 0 ? g_token : stage_token(); }
// g_host_is_amd stands in for the host CPU vendor, so the AMD refusal is
// tested on every host and the other 1 GB cases run on AMD hosts too;
// g_rmid_forced likewise for kernel.shm_rmid_forced (1 set, 0 not, -1 it
// cannot be read). g_rmid_fails makes the removal of a System V segment fail
// with EPERM, as for another user's.
static bool g_host_is_amd = false;
static bool test_host_is_amd(void) { return g_host_is_amd; }
#define BWA_SHM_HOST_IS_AMD test_host_is_amd
static int g_rmid_forced = 0;
static int test_rmid_forced(void) {
    if (g_rmid_forced < 0) errno = EACCES;
    return g_rmid_forced;
}
#define BWA_SHM_RMID_FORCED test_rmid_forced
static bool g_rmid_fails = false;
static int test_rmid(int id) {
    if (!g_rmid_fails) return shmctl(id, IPC_RMID, NULL);
    errno = EPERM;
    return -1;
}
#define BWA_SHM_RMID(id) test_rmid(id)
#define BWA_SHM_HUGE_FLAGS g_huge_flags
#define BWA_SHM_HUGE_PAGE g_huge_page
#define BWA_SHM_STAGE_TOKEN test_token
#define BWA_SHM_STAGE_WAIT_S g_wait_s
#define BWA_SHM_STATVFS fake_statvfs
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
struct Names {
    bwa_tests::TempDir dir{"ss"};
    std::string segment;                 // the first base's segment
    std::vector<std::string> segments;
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
    }
    ~Names() {
        for (const std::string &seg : segments) shm_unlink(seg.c_str());
        shm_unlink(private_control);
        sem_unlink(private_lock);
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

// The stager lock exists.
bool lock_exists() {
    sem_t *s = sem_open(private_lock, 0);
    if (s == SEM_FAILED) return false;
    sem_close(s);
    return true;
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

// Release what bwa_shm_attach returned, as the aligner does.
void detach(uint8_t *m, size_t map_len) { bwa_shm_detach(m, map_len); }

// Attach succeeds and its bytes equal bwa_shm_pack_from_disk's; with
// `want_map_len`, the mapping is that long.
void check_attach_matches(const std::string &prefix, size_t want_map_len = 0) {
    size_t len = 0, map_len = 0;
    uint8_t *m = bwa_shm_attach(prefix.c_str(), &len, &map_len);
    CHECK(m != NULL);
    CHECK(map_len >= len);
    CHECK(want_map_len == 0 || map_len == want_map_len);
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
    for (uint8_t *v : {a.data(), b.data()}) {
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
    detach(m, map_len);
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
// When set, a stage child's stderr goes to this file (see stage_stderr).
std::string g_child_stderr;
// Every child this test forked that may create a System V segment: the
// teardown of a 1 GB case removes only segments one of them created.
std::vector<pid_t> g_children;

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
pid_t spawn_stage(const std::string &prefix, HookMode mode, bool huge = false) {
    fflush(NULL);
    const pid_t pid = fork();
    CHECK(pid >= 0);
    if (pid == 0) {
        for (int fd : g_child_close) close(fd);
        alarm(20);
        if (!g_child_stderr.empty() && freopen(g_child_stderr.c_str(), "w", stderr) == NULL)
            _exit(94);
        g_mode = mode;
        bwa_shm_stage_map_hook = stage_hook;
        const int rc = test_private_shm_stage(prefix.c_str(), false, -1, 1, huge);
        fflush(NULL);
        _exit(rc == 0 ? 0 : 1);
    }
    g_children.push_back(pid);
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

// The host is Linux, where kernel.shm_rmid_forced exists to be read.
bool is_linux() {
#ifdef __linux__
    return true;
#else
    return false;
#endif
}

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
    explicit PausedStager(const std::string &prefix, bool huge = false) {
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

} // namespace

// ---- 1 GB staging (`shm --huge-pages`): the index in a System V segment ----

// One case's 1 GB stages. Each stage takes its token from fresh(), so its key
// is known in this process: if a CHECK throws, destruction removes what the
// case staged, the System V segments included, and puts the test seams back.
// A segment on one of those keys is removed only if one of this test's own
// children created it (g_children): the key alone is no evidence, since
// another program's segment may hold it. Declared after the case's Names, so
// it runs while the registry still exists.
struct HugeCase {
    const uint64_t base = stage_token();
    std::vector<uint64_t> tokens;
    std::vector<int> own_ids;   // segments a child made that IPC_STAT cannot read
    uint64_t fresh() {
        g_token = base + 2 * (tokens.size() + 1);
        tokens.push_back(g_token);
        return g_token;
    }
    ~HugeCase() {
        g_rmid_fails = g_host_is_amd = false;
        g_rmid_forced = 0;
        g_fake_free_bytes = -1;
        g_fake_frsize = 1;
        g_destroy_wait_ms = 5000;
        g_child_stderr.clear();
        unsetenv("BWAMEM3_SHM_HUGE_PAGES_AMD");
        test_private_shm_destroy();
        for (const uint64_t t : tokens) {
            const int id = shmget(huge_key(t), 0, 0);
            struct shmid_ds ds;
            if (id >= 0 && shmctl(id, IPC_STAT, &ds) == 0 &&
                std::find(g_children.begin(), g_children.end(), ds.shm_cpid) != g_children.end())
                shmctl(id, IPC_RMID, NULL);
        }
        for (const int id : own_ids) shmctl(id, IPC_RMID, NULL);
        g_token = 0;
        g_huge_flags = 0;
        g_huge_page = 65536;
    }
};

// A System V segment made by this test process: "another program's".
struct ForeignSegment {
    key_t key;
    int id;
    ForeignSegment(uint64_t token, size_t bytes) : key(huge_key(token)) {
        id = shmget(key, bytes, IPC_CREAT | IPC_EXCL | 0644);
        CHECK(id >= 0);
    }
    bool exists() const { return shmget(key, 0, 0) == id; }
    ~ForeignSegment() { shmctl(id, IPC_RMID, NULL); }
};

bool key_absent(uint64_t token) {
    return shmget(huge_key(token), 0, 0) < 0 && errno == ENOENT;
}

// The PID of a process that has exited.
uint64_t dead_pid() {
    const pid_t pid = fork();
    CHECK(pid >= 0);
    if (pid == 0) _exit(0);
    wait_status(pid);
    return (uint64_t)pid;
}

// A System V segment of `bytes` bytes with `token`'s key and permissions
// `mode`, made by a child that has exited: what a stager leaves. Returns the
// child's PID.
uint64_t dead_stagers_segment(uint64_t token, size_t bytes, int mode = 0600) {
    fflush(NULL);
    const pid_t stager = fork();
    CHECK(stager >= 0);
    if (stager == 0)
        _exit(shmget(huge_key(token), bytes, IPC_CREAT | IPC_EXCL | mode) >= 0 ? 0 : 1);
    g_children.push_back(stager);
    CHECK(exited(wait_status(stager), 0) && !key_absent(token));
    return (uint64_t)stager;
}

// Set reserved[`i`] in the header of POSIX segment `name`.
void set_reserved(const std::string &name, int i, uint64_t value) {
    const int fd = shm_open(name.c_str(), O_RDWR, 0);
    CHECK(fd >= 0);
    void *m = mmap(NULL, sizeof(bwa_shm_header_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
    CHECK(m != MAP_FAILED);
    ((bwa_shm_header_t *)m)->reserved[i] = value;
    munmap(m, sizeof(bwa_shm_header_t));
}

// Set the version in the header of POSIX segment `name`.
void set_version(const std::string &name, uint32_t version) {
    const int fd = shm_open(name.c_str(), O_RDWR, 0);
    CHECK(fd >= 0);
    void *m = mmap(NULL, sizeof(bwa_shm_header_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
    CHECK(m != MAP_FAILED);
    ((bwa_shm_header_t *)m)->version = version;
    munmap(m, sizeof(bwa_shm_header_t));
}

// Remove registry entry `name` by hand (under the private lock).
void unregister_entry(const char *name) {
    sem_t *lock = ctl_lock_acquire();
    CHECK(lock != NULL);
    int fd = -1;
    void *m = ctl_open_rw(&fd);
    CHECK(m != NULL);
    ctl_remove(m, name);
    ctl_close(m, fd);
    ctl_lock_release(lock);
}

// This process's `kind` ("ipc", "pid") namespace, read independently of the
// implementation: the number in the link text "<kind>:[N]". 0 where there is
// none (not Linux).
uint64_t namespace_from_link(const char *kind) {
    char buf[64];
    const ssize_t n = readlink((std::string("/proc/self/ns/") + kind).c_str(), buf,
                               sizeof buf - 1);
    if (n <= 0) return 0;
    buf[n] = '\0';
    unsigned long long ns = 0;
    CHECK(sscanf(buf, (std::string(kind) + ":[%llu]").c_str(), &ns) == 1 && ns != 0);
    return (uint64_t)ns;
}

// A 1 GB stub made by hand: `magic`, stager `pid`, `token`, and a System V
// segment of `bytes` bytes with no id and no IPC namespace recorded.
void make_stub(const std::string &name, uint64_t magic, uint64_t pid, uint64_t token,
               uint64_t bytes) {
    make_segment(name, sizeof(bwa_shm_header_t), magic, pid, token, /*declared=*/4096);
    set_reserved(name, 2, bytes);
}

std::string slurp(const std::string &path) {
    std::string out;
    FILE *f = fopen(path.c_str(), "rb");
    CHECK(f != NULL);
    char buf[4096];
    size_t got;
    while ((got = fread(buf, 1, sizeof buf, f)) > 0) out.append(buf, got);
    fclose(f);
    return out;
}

bool has(const std::string &text, const char *what) {
    return text.find(what) != std::string::npos;
}

// Stage `prefix` in a child that must exit with `code`; returns its stderr.
std::string stage_stderr(const Names &n, const std::string &prefix, bool huge, int code) {
    g_child_stderr = n.dir.file("stderr");
    const int status = wait_status(spawn_stage(prefix, HOOK_NONE, huge));
    const std::string path = g_child_stderr;
    g_child_stderr.clear();
    CHECK(exited(status, code));
    return slurp(path);
}

// The stage was killed mid-pack and left what a 1 GB stage must leave: an
// incomplete stub that names the segment its stager created. A stage that died
// of another signal, or was not a 1 GB stage, does not count.
bool killed_huge_stage(int status, const std::string &segment) {
    if (!WIFSIGNALED(status) || WTERMSIG(status) != SIGKILL) return false;
    const bwa_shm_header_t stub = read_header(segment);
    return stub.magic == 0 && is_huge_stub(&stub) && huge_segment(&stub) >= 0;
}

// h1. A 1 GB stage leaves a stub and the index in the System V segment the
//     stub names, of whole pages; attach follows it; a restage changes nothing
//     and warns about the other page kind; a stub whose segment was removed is
//     reclaimed; -d removes both.
void check_huge_stage() {
    uint64_t total = 0;
    {
        Names n({"size"});
        uint8_t *disk = NULL;
        CHECK(bwa_shm_pack_from_disk(n.prefix("size").c_str(), &disk, &total) == 0);
        free(disk);
    }
    CHECK(total > 16384 && total < 32768);
    // One page, two pages, and an index that fills its page exactly.
    for (const unsigned long long page : {65536ull, 16384ull, (unsigned long long)total}) {
        Names n({"huge", "plain"});
        HugeCase c;
        g_huge_page = page;
        const uint64_t want = (total + page - 1) / page * page;
        const std::string prefix = n.prefix("huge");
        c.fresh();
        CHECK(!has(stage_stderr(n, prefix, /*huge=*/true, 0), "[W::"));
        const bwa_shm_header_t stub = read_header(n.segment);
        CHECK(stub.magic == BWA_SHM_HUGE_MAGIC && is_huge_stub(&stub));
        CHECK(stub.total_size == total && stub.reserved[2] == want);
        // The stager's IPC namespace in the low half, its PID namespace in the high.
        CHECK(stub.reserved[4] ==
              (namespace_from_link("ipc") | namespace_from_link("pid") << 32));
        const int id = huge_segment(&stub);
        CHECK(id >= 0 && stub.reserved[3] == (uint64_t)id + 1);
        struct shmid_ds ds;
        CHECK(shmctl(id, IPC_STAT, &ds) == 0 && (uint64_t)ds.shm_segsz == want);
        CHECK((uint64_t)ds.shm_cpid == stub.reserved[0] && (ds.shm_perm.mode & 0777) == 0644);
        check_attach_matches(prefix, (size_t)want);
        // The attachment was given back whole.
        CHECK(shmctl(id, IPC_STAT, &ds) == 0 && ds.shm_nattch == 0);
        CHECK(list_output() == "huge\t" + std::to_string(total) + "\n");
        // Restaging keeps the same stub and segment, and without --huge-pages
        // says which kind of pages the index is on.
        for (const bool huge : {true, false}) {
            c.fresh();
            const std::string err = stage_stderr(n, prefix, huge, 0);
            CHECK(has(err, "is already in shm") && !has(err, "[W::") && !has(err, "[E::"));
            CHECK(has(err, "is already in shm (on 1 GB pages)") == !huge);
            const bwa_shm_header_t again = read_header(n.segment);
            CHECK(std::memcmp(&again, &stub, sizeof stub) == 0 && huge_segment(&again) == id);
        }
        // The reverse: an index staged in POSIX shm stays there, and the 1 GB
        // stage that asked for other pages fails.
        const std::string plain = n.prefix("plain"), plain_seg = n.segments.at(1);
        CHECK(!has(stage_stderr(n, plain, /*huge=*/false, 0), "[W::"));
        const bwa_shm_header_t before = read_header(plain_seg);
        c.fresh();
        CHECK(has(stage_stderr(n, plain, /*huge=*/true, 1),
                  "is already in shm, and not on 1 GB pages"));
        CHECK(key_absent(g_token));
        const bwa_shm_header_t after = read_header(plain_seg);
        CHECK(before.magic == BWA_SHM_MAGIC && std::memcmp(&before, &after, sizeof before) == 0);
        // The segment removed behind the stub (ipcrm, or the host): attach
        // falls back, the listing says incomplete, and a restage reclaims it.
        CHECK(shmctl(id, IPC_RMID, NULL) == 0);
        CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
        CHECK(has(list_output(), ("huge\t" + std::to_string(total) + "\tincomplete\n").c_str()));
        const uint64_t token = c.fresh();
        CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 0),
                  "reclaiming 'huge': its 1 GB segment"));
        const bwa_shm_header_t restaged = read_header(n.segment);
        CHECK(restaged.magic == BWA_SHM_HUGE_MAGIC && restaged.reserved[1] == token);
        CHECK(huge_segment(&restaged) >= 0);
        check_attach_matches(prefix, (size_t)want);
        CHECK(test_private_shm_destroy() == 0);
        CHECK(!segment_exists(n.segment) && !segment_exists(plain_seg));
        CHECK(key_absent(token));
        CHECK(shm_open(private_control, O_RDONLY, 0) < 0 && errno == ENOENT);
    }
    // The key is never IPC_PRIVATE, which would make a segment nothing can find.
    CHECK(IPC_PRIVATE == 0 && huge_key(1) == 1 && huge_key(2) == 1 && huge_key(4) == 2);
    std::printf("ok h1 1 GB stage: stub and System V segment; attach, restage, destroy\n");
}

// h2. A 1 GB stage killed mid-pack leaves its segment named by a registered
//     stub, so a restage, or `shm -d`, removes it.
void check_huge_interrupted() {
    for (const bool restage : {true, false}) {
        Names n({"hkill"});
        HugeCase c;
        const std::string prefix = n.prefix("hkill");
        // An ordinary stage killed at the same point is not a 1 GB leftover.
        CHECK(!killed_huge_stage(wait_status(spawn_stage(prefix, HOOK_SIGKILL)), n.segment));
        CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
        CHECK(test_private_shm_destroy() == 0);
        const uint64_t token = c.fresh();
        CHECK(killed_huge_stage(wait_status(spawn_stage(prefix, HOOK_SIGKILL, true)),
                                n.segment));
        CHECK(read_header(n.segment).reserved[1] == token && !key_absent(token));
        CHECK(lock_free());
        CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
        CHECK(list_output().find("hkill\t") == 0);
        CHECK(has(list_output(), "\tincomplete\n"));
        if (restage) {
            c.fresh();
            CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 0), "left by an interrupted stage"));
            CHECK(key_absent(token));
            check_attach_matches(prefix);
            CHECK(entries() == std::vector<std::string>{"hkill"});
        }
        CHECK(test_private_shm_destroy() == 0);
        CHECK(key_absent(token) && key_absent(g_token) && !segment_exists(n.segment));
    }
    std::printf("ok h2 interrupted 1 GB stage: segment removed by a restage and by destroy\n");
}

// h3. A System V segment the stage did not create is never removed: not when
//     it already holds the stage's key, and not when a stub names it. One that
//     cannot be seen from this IPC namespace is not taken for gone.
void check_huge_foreign() {
    Names n({"hforeign"});
    HugeCase c;
    const std::string prefix = n.prefix("hforeign");
    const uint64_t token = c.fresh();
    ForeignSegment foreign(token, 65536);
    // The stage's key is taken.
    CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 1), "another segment already has"));
    CHECK(foreign.exists());
    CHECK(entries().empty() && !segment_exists(n.segment) && lock_free());
    // Complete stubs with the foreign segment's key that their stage did not
    // create: another creator; a recorded id that is not the segment's; the
    // right creator but another size. Each is unusable, not attached, not
    // replaced, and destroy removes only the stub.
    const uint64_t self = (uint64_t)getpid();
    const struct { uint64_t pid, id_plus_1, bytes; } stubs[] = {
        {dead_pid(), 0, 65536}, {self, (uint64_t)foreign.id + 2, 65536}, {self, 0, 131072}};
    for (const auto &s : stubs) {
        make_stub(n.segment, BWA_SHM_HUGE_MAGIC, s.pid, token, s.bytes);
        set_reserved(n.segment, 3, s.id_plus_1);
        register_entry("hforeign", 4096);
        CHECK(has(list_output(), "\tunusable\n"));
        CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
        CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 1), "it is another program's"));
        CHECK(read_header(n.segment).magic == BWA_SHM_HUGE_MAGIC);   // not replaced
        CHECK(test_private_shm_destroy() == 0);
        CHECK(foreign.exists() && !segment_exists(n.segment));
    }
    // A stager killed after creating its segment and before recording the id:
    // the segment is the one the stub's PID created, and destroy removes it.
    const uint64_t token2 = c.fresh();
    const uint64_t stager = dead_stagers_segment(token2, 65536);
    make_stub(n.segment, 0, stager, token2, 65536);
    register_entry("hforeign", 4096);
    CHECK(has(list_output(), "\tincomplete\n"));
    CHECK(test_private_shm_destroy() == 0);
    CHECK(key_absent(token2) && foreign.exists() && !segment_exists(n.segment));
    // The same, but the stager was in another PID namespace: the creator PID
    // the kernel reports here cannot be compared with the stub's, so the stub
    // and the segment are kept (destroy fails) until seen from that namespace.
    const uint64_t token4 = c.fresh();
    const uint64_t stager4 = dead_stagers_segment(token4, 65536);
    make_stub(n.segment, 0, stager4, token4, 65536);
    set_reserved(n.segment, 4, ipc_namespace() | (pid_namespace() + 1) << 32);
    register_entry("hforeign", 4096);
    CHECK(test_private_shm_destroy() != 0);
    CHECK(segment_exists(n.segment) && !key_absent(token4));
    CHECK(entries() == std::vector<std::string>{"hforeign"});
    set_reserved(n.segment, 4, stub_namespaces());
    CHECK(test_private_shm_destroy() == 0);
    CHECK(key_absent(token4) && !segment_exists(n.segment));
    // A segment on the stub's key that may not be read (IPC_STAT denied) is
    // not taken for another program's: that does not show whose it is, so the
    // stub stays to name it. (Root reads every segment.)
    if (geteuid() == 0) {
        std::printf("SKIP: h3 segment that may not be inspected (running as root)\n");
    } else {
        const uint64_t token5 = c.fresh();
        const uint64_t stager5 = dead_stagers_segment(token5, 65536, /*mode=*/0);
        c.own_ids.push_back(shmget(huge_key(token5), 0, 0));   // unreadable, so by id
        make_stub(n.segment, 0, stager5, token5, 65536);
        register_entry("hforeign", 4096);
        CHECK(test_private_shm_destroy() != 0);
        CHECK(segment_exists(n.segment) && !key_absent(token5));
        CHECK(entries() == std::vector<std::string>{"hforeign"});
        const int id5 = shmget(huge_key(token5), 0, 0);
        CHECK(id5 >= 0 && shmctl(id5, IPC_RMID, NULL) == 0);   // its owner may remove it
        CHECK(test_private_shm_destroy() == 0);
        CHECK(key_absent(token5) && !segment_exists(n.segment));
    }
    // A complete stub staged in another IPC namespace: its segment's key means
    // nothing here, so it is unusable and neither a stage nor destroy touches
    // the stub or the registry. Seen from its own namespace it is complete.
    const uint64_t token3 = c.fresh();
    const uint64_t stager3 = dead_stagers_segment(token3, 65536);
    make_stub(n.segment, BWA_SHM_HUGE_MAGIC, stager3, token3, 65536);
    const int id3 = shmget(huge_key(token3), 0, 0);
    set_reserved(n.segment, 3, (uint64_t)id3 + 1);
    // The segment holds a complete index header with the stage's token, as a
    // completed stage leaves it.
    bwa_shm_header_t *index3 = (bwa_shm_header_t *)shmat(id3, NULL, 0);
    CHECK(index3 != (void *)-1);
    index3->version = BWA_SHM_VERSION;
    index3->total_size = 4096;
    index3->reserved[1] = token3;
    index3->magic = BWA_SHM_MAGIC;
    CHECK(shmdt(index3) == 0);
    set_reserved(n.segment, 4, ipc_namespace() + 1);
    register_entry("hforeign", 4096);
    CHECK(has(list_output(), "\tunusable\n"));
    CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
    // (Where a process cannot read its namespace, the message says that.)
    CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 1),
              ipc_namespace() != 0 ? "it is in another IPC namespace" : "its IPC namespace"));
    CHECK(test_private_shm_destroy() != 0);
    CHECK(segment_exists(n.segment) && !key_absent(token3));
    CHECK(entries() == std::vector<std::string>{"hforeign"});
    set_reserved(n.segment, 4, ipc_namespace());
    CHECK(list_output() == "hforeign\t4096\n");
    CHECK(test_private_shm_destroy() == 0);
    CHECK(key_absent(token3) && foreign.exists() && !segment_exists(n.segment));
    std::printf("ok h3 only the System V segment a stub's own stager created is removed\n");
}

// h4. A 1 GB stage that fails, loses its stub to `shm -d` or to another
//     stager, or exits mid-pack, removes the segment it created.
void check_huge_rollback() {
    for (const HookMode mode : {HOOK_FIFO_INPUT, HOOK_DESTROY, HOOK_FATAL_INPUT}) {
        Names n({"hfail"});
        HugeCase c;
        const std::string prefix = n.prefix("hfail", /*copy=*/true);
        g_break_path = prefix + (mode == HOOK_FATAL_INPUT ? ".pac" : ".bwt.2bit.64");
        const uint64_t token = c.fresh();
        CHECK(exited(wait_status(spawn_stage(prefix, mode, /*huge=*/true)), 1));
        CHECK(key_absent(token));
        if (mode != HOOK_FATAL_INPUT) {
            CHECK(entries().empty() && !segment_exists(n.segment));
            continue;
        }
        // The pack exited the process, so nothing rolled the stub back: it
        // stays, incomplete, with no segment behind it.
        const bwa_shm_header_t stub = read_header(n.segment);
        CHECK(stub.magic == 0 && is_huge_stub(&stub) && stub.reserved[1] == token);
        CHECK(has(list_output(), "\tincomplete\n"));
        CHECK(test_private_shm_destroy() == 0);
        CHECK(!segment_exists(n.segment));
    }
    {
        // The name was taken over mid-stage (as in 7): the replacement stays,
        // and the segment this stage created still goes.
        Names n({"htaken"});
        HugeCase c;
        const std::string prefix = n.prefix("htaken");
        const uint64_t token = c.fresh();
        PausedStager a(prefix, /*huge=*/true);
        CHECK(!key_absent(token));
        CHECK(shm_unlink(n.segment.c_str()) == 0);
        make_segment(n.segment, sizeof(bwa_shm_header_t), 0, (uint64_t)getpid(), 2);
        CHECK(exited(a.finish(), 1));
        CHECK(key_absent(token));
        const bwa_shm_header_t h = read_header(n.segment);
        CHECK(h.magic == 0 && h.reserved[1] == 2);
        CHECK(entries() == std::vector<std::string>{"htaken"});
    }
#ifdef __linux__
    {
        // The segment removed by hand mid-stage (ipcrm): it lives on, marked
        // for removal, only while the stager is attached, so the stage fails
        // rather than publish an index that is about to vanish.
        Names n({"hipcrm"});
        HugeCase c;
        const std::string prefix = n.prefix("hipcrm");
        const uint64_t token = c.fresh();
        PausedStager a(prefix, /*huge=*/true);
        const int id = shmget(huge_key(token), 0, 0);
        CHECK(id >= 0 && shmctl(id, IPC_RMID, NULL) == 0);
        CHECK(exited(a.finish(), 1));
        CHECK(entries().empty() && !segment_exists(n.segment));
    }
#endif
    std::printf("ok h4 failed, destroyed, taken-over or exited 1 GB stage leaves no segment\n");
}

// h5. A host that cannot take a 1 GB stage refuses one and stages nothing: an
//     AMD x86 host unless BWAMEM3_SHM_HUGE_PAGES_AMD=1, one whose kernel
//     removes a segment when its stager exits or whose setting for that cannot
//     be read, and one that is not Linux. An ordinary stage is not affected.
void check_huge_refused() {
    Names n({"hno"});
    HugeCase c;
    const std::string prefix = n.prefix("hno");
    const auto refused = [&](const char *why) {
        c.fresh();
        CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 1), why));
        CHECK(entries().empty() && !segment_exists(n.segment) && key_absent(g_token));
    };
    CHECK(unsetenv("BWAMEM3_SHM_HUGE_PAGES_AMD") == 0);
    g_host_is_amd = true;
    for (const char *value : {(const char *)NULL, "0", "", "yes"}) {
        if (value != NULL) CHECK(setenv("BWAMEM3_SHM_HUGE_PAGES_AMD", value, 1) == 0);
        refused("--huge-pages is refused on AMD x86 hosts");
    }
    CHECK(setenv("BWAMEM3_SHM_HUGE_PAGES_AMD", "1", 1) == 0);
    c.fresh();
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, /*huge=*/true)), 0));
    CHECK(unsetenv("BWAMEM3_SHM_HUGE_PAGES_AMD") == 0);
    CHECK(read_header(n.segment).magic == BWA_SHM_HUGE_MAGIC && !key_absent(g_token));
    CHECK(test_private_shm_destroy() == 0);
    CHECK(key_absent(g_token));
    // An ordinary stage is allowed, and a 1 GB stage of that index is still
    // refused, not reported as already staged.
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE)), 0));
    c.fresh();
    CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 1), "refused on AMD x86 hosts"));
    CHECK(read_header(n.segment).magic == BWA_SHM_MAGIC && key_absent(g_token));
    CHECK(test_private_shm_destroy() == 0);
    g_host_is_amd = false;
    g_rmid_forced = 1;
    refused("kernel.shm_rmid_forced=1");
    g_rmid_forced = -1;   // unknown is not taken for 0
    refused("--huge-pages cannot read kernel.shm_rmid_forced");
    g_rmid_forced = 0;
    g_huge_flags = -1;
    refused("--huge-pages needs Linux");
    g_huge_flags = 0;
    std::printf("ok h5 1 GB stage refused on AMD x86 (unless overridden), with "
                "shm_rmid_forced, and off Linux\n");
}

// h6. A 1 GB segment that cannot be removed keeps its stub, its entry and the
//     registry, and `shm -d` and a reclaiming stage fail, so a later attempt
//     finds it; the entries of the segments that were removed go.
void check_huge_removal_fails() {
    Names n({"hkeep", "hgo", "hstuck"});
    HugeCase c;
    const std::string keep = n.prefix("hkeep"), go = n.prefix("hgo");
    const uint64_t token = c.fresh();
    CHECK(exited(wait_status(spawn_stage(keep, HOOK_NONE, /*huge=*/true)), 0));
    CHECK(exited(wait_status(spawn_stage(go, HOOK_NONE)), 0));
    g_rmid_fails = true;
    CHECK(test_private_shm_destroy() != 0);
    CHECK(entries() == std::vector<std::string>{"hkeep"});
    CHECK(segment_exists(n.segments[0]) && !segment_exists(n.segments[1]));
    CHECK(!key_absent(token) && lock_exists() && lock_free());
    check_attach_matches(keep);
    g_rmid_fails = false;
    CHECK(test_private_shm_destroy() == 0);
    CHECK(key_absent(token) && !segment_exists(n.segments[0]));
    // The reclaim of a dead 1 GB stage that cannot remove its segment.
    const uint64_t dead = c.fresh();
    CHECK(killed_huge_stage(wait_status(spawn_stage(keep, HOOK_SIGKILL, true)), n.segment));
    g_rmid_fails = true;
    c.fresh();
    CHECK(has(stage_stderr(n, keep, /*huge=*/true, 1), "cannot remove the 1 GB segment"));
    CHECK(!key_absent(dead) && read_header(n.segment).reserved[1] == dead);
    CHECK(entries() == std::vector<std::string>{"hkeep"} && lock_free());
    g_rmid_fails = false;
    c.fresh();
    CHECK(exited(wait_status(spawn_stage(keep, HOOK_NONE, /*huge=*/true)), 0));
    CHECK(key_absent(dead));
    check_attach_matches(keep);
    CHECK(test_private_shm_destroy() == 0);
    // A stage that fails after creating its segment, and cannot remove it,
    // rolls back neither the stub that names it nor the entry.
    const std::string stuck = n.prefix("hstuck", /*copy=*/true);
    g_break_path = stuck + ".bwt.2bit.64";
    const uint64_t stuck_token = c.fresh();
    g_rmid_fails = true;
    CHECK(exited(wait_status(spawn_stage(stuck, HOOK_FIFO_INPUT, /*huge=*/true)), 1));
    g_rmid_fails = false;
    CHECK(!key_absent(stuck_token) && segment_exists(n.segments[2]));
    CHECK(read_header(n.segments[2]).reserved[1] == stuck_token);
    CHECK(entries() == std::vector<std::string>{"hstuck"});
    CHECK(test_private_shm_destroy() == 0);
    CHECK(key_absent(stuck_token) && !segment_exists(n.segments[2]));
    std::printf("ok h6 an unremovable 1 GB segment keeps its stub and entry\n");
}

// h7. Attach maps the segment only if it has the size the stub records and
//     holds the index that the stub's stage packed; `shm -l` and a stage judge
//     it the same way.
void check_huge_attach_validation() {
    Names n({"hval"});
    HugeCase c;
    g_huge_page = 16384;   // two pages, so the segment can also be too large
    const std::string prefix = n.prefix("hval");
    const uint64_t token = c.fresh();
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, /*huge=*/true)), 0));
    const bwa_shm_header_t stub = read_header(n.segment);
    const int id = huge_segment(&stub);
    CHECK(id >= 0);
    CHECK(stub.reserved[2] == 2 * g_huge_page);
    // A stub that records more bytes than the segment has (mapping that many
    // would run past it), or fewer.
    // Listed and judged by a stage as attach judges it: unusable, and not
    // replaced, rather than reported as already staged.
    for (const uint64_t bytes : {stub.reserved[2] + g_huge_page, stub.reserved[2] - g_huge_page}) {
        set_reserved(n.segment, 2, bytes);
        CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
        CHECK(has(list_output(), "hval\t") && has(list_output(), "\tunusable\n"));
        c.fresh();
        const std::string err = stage_stderr(n, prefix, /*huge=*/true, 1);
        CHECK(has(err, "it is not the size its stub records") && !has(err, "already in shm"));
        CHECK(has(err, "bwa-mem3 shm -d removes it and its segment"));
        CHECK(read_header(n.segment).reserved[2] == bytes && huge_segment(&stub) == id);
    }
    set_reserved(n.segment, 2, stub.reserved[2]);
    check_attach_matches(prefix);
    // A segment whose index carries another stage's token, or is not complete.
    bwa_shm_header_t *h = (bwa_shm_header_t *)shmat(id, NULL, 0);
    CHECK(h != (void *)-1 && h->reserved[1] == token && h->magic == BWA_SHM_MAGIC);
    h->reserved[1] = token + 2;
    uint8_t *m = bwa_shm_attach(prefix.c_str(), NULL);
    const bool listed_unusable = has(list_output(), "\tunusable\n");
    h->reserved[1] = token;
    CHECK(m == NULL && listed_unusable);
    h->magic = 0;
    m = bwa_shm_attach(prefix.c_str(), NULL);
    const bool incomplete_unusable = has(list_output(), "\tunusable\n");
    h->magic = BWA_SHM_MAGIC;
    CHECK(incomplete_unusable);
    // An index header of another version, or declaring more bytes than the
    // segment holds.
    for (const int field : {0, 1}) {
        const uint32_t version = h->version;
        const uint64_t total = h->total_size;
        if (field == 0) h->version = BWA_SHM_VERSION + 1;
        else h->total_size = stub.reserved[2] + 1;
        m = bwa_shm_attach(prefix.c_str(), NULL);
        const bool unusable = has(list_output(), "\tunusable\n");
        c.fresh();
        const std::string err = stage_stderr(n, prefix, /*huge=*/true, 1);
        h->version = version;
        h->total_size = total;
        CHECK(m == NULL && unusable);
        CHECK(has(err, "another version or a size past the segment"));
    }
    CHECK(shmdt(h) == 0);
    CHECK(m == NULL);
    check_attach_matches(prefix);
    // A stub of another version is not attached and not replaced, and the
    // message names the version; `shm -d` still removes it and its segment.
    set_version(n.segment, BWA_SHM_VERSION + 1);
    CHECK(bwa_shm_attach(prefix.c_str(), NULL) == NULL);
    CHECK(has(list_output(), "\tunusable\n"));
    c.fresh();
    CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 1),
              ("of version " + std::to_string(BWA_SHM_VERSION + 1)).c_str()));
    CHECK(huge_segment(&stub) == id);
    set_version(n.segment, BWA_SHM_VERSION);
    // A complete stub that lost its registry entry is not replaced either.
    unregister_entry("hval");
    c.fresh();
    CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 1), "exists without a registry entry"));
    CHECK(huge_segment(&stub) == id && entries().empty());
    register_entry("hval", (int64_t)stub.total_size);
    check_attach_matches(prefix);
    CHECK(test_private_shm_destroy() == 0);
    CHECK(key_absent(token));
    std::printf("ok h7 attach checks the 1 GB segment's size and token\n");
}

// h8. /dev/shm needs room for a 1 GB stage's stub only, and a 1 GB stage in
//     progress reserves nothing there for the stages beside it.
void check_huge_preflight() {
    Names n({"hfull", "hstub", "hposix"});
    HugeCase c;
    const std::string prefix = n.prefix("hfull"), posix = n.prefix("hposix");
    g_fake_free_bytes = (long long)sizeof(bwa_shm_header_t) - 1;
    c.fresh();
    CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 1), "no room for the 64-byte stub"));
    CHECK(entries().empty() && !segment_exists(n.segment) && key_absent(g_token));
    CHECK(lock_free());
    g_fake_free_bytes = (long long)sizeof(bwa_shm_header_t);
    c.fresh();
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, /*huge=*/true)), 0));
    CHECK(test_private_shm_destroy() == 0);
    uint8_t *disk = NULL;
    uint64_t need = 0;
    CHECK(bwa_shm_pack_from_disk(posix.c_str(), &disk, &need) == 0);
    free(disk);
    make_stub(n.segments[1], 0, (uint64_t)getpid(), c.fresh(), 65536);   // a live stager's
    register_entry("hstub", (int64_t)need);
    g_fake_free_bytes = (long long)need;
    CHECK(exited(wait_status(spawn_stage(posix, HOOK_NONE)), 0));
    g_fake_free_bytes = -1;
    CHECK(test_private_shm_destroy() == 0);
    std::printf("ok h8 /dev/shm preflight counts only the stub of a 1 GB stage\n");
}

// h8b. Every segment takes whole tmpfs pages: a stub is refused when the only
//      free page is the one another stage in progress still needs for the
//      partial last page of its index.
void check_huge_preflight_last_page() {
    Names n({"hlast", "hpart", "hpart2"});
    HugeCase c;
    const std::string prefix = n.prefix("hlast");
    const unsigned long unit = 4096;
    // A live ordinary stage, 64 bytes into a page it has not yet written.
    const auto in_progress = [&](int i, const char *name) {
        make_segment(n.segments[i], sizeof(bwa_shm_header_t), 0, (uint64_t)getpid(), 9);
        const int fd = shm_open(n.segments[i].c_str(), O_RDONLY, 0);
        CHECK(fd >= 0);
        struct stat st;
        CHECK(fstat(fd, &st) == 0);
        close(fd);
        const uint64_t allocated = st.st_blocks > 0 ? (uint64_t)st.st_blocks * 512u : 0;
        CHECK(allocated % unit == 0);
        register_entry(name, (int64_t)(allocated + 64));
    };
    in_progress(1, "hpart");
    g_fake_frsize = unit;
    g_fake_free_bytes = (long long)unit;
    c.fresh();
    CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 1), "no room for the 64-byte stub"));
    CHECK(entries() == std::vector<std::string>{"hpart"});
    CHECK(!segment_exists(n.segments[0]) && key_absent(g_token));
    // Two such stages need a page each: with two free pages none is left for
    // the stub, though their unwritten bytes alone would leave nearly two.
    in_progress(2, "hpart2");
    g_fake_free_bytes = 2 * (long long)unit;
    c.fresh();
    CHECK(has(stage_stderr(n, prefix, /*huge=*/true, 1), "no room for the 64-byte stub"));
    CHECK(!segment_exists(n.segments[0]) && key_absent(g_token));
    g_fake_free_bytes = 3 * (long long)unit;
    c.fresh();
    CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, /*huge=*/true)), 0));
    CHECK(read_header(n.segments[0]).magic == BWA_SHM_HUGE_MAGIC);
    g_fake_free_bytes = -1;
    CHECK(test_private_shm_destroy() == 0);
    std::printf("ok h8b /dev/shm preflight counts whole pages\n");
}

// h10. A stage whose System V segment could not be created rolls back only its
//      own stub and entry. Here one process stages two indexes of one size
//      that draw the same key: the second stage's shmget fails, and the first
//      index's segment, which the second stub's PID and size also match, stays.
void check_huge_rollback_same_key() {
    Names n({"hdup1", "hdup2"});
    HugeCase c;
    const std::string first = n.prefix("hdup1"), second = n.prefix("hdup2");
    const uint64_t token = c.fresh();
    const std::string path = n.dir.file("stderr");
    fflush(NULL);
    const pid_t pid = fork();
    CHECK(pid >= 0);
    if (pid == 0) {
        alarm(20);
        if (freopen(path.c_str(), "w", stderr) == NULL) _exit(94);
        const int a = test_private_shm_stage(first.c_str(), false, -1, 1, true);
        const int b = test_private_shm_stage(second.c_str(), false, -1, 1, true);
        fflush(NULL);
        _exit(a == 0 && b != 0 ? 0 : 1);
    }
    g_children.push_back(pid);
    CHECK(exited(wait_status(pid), 0));
    CHECK(has(slurp(path), "another segment already has this stage's System V key"));
    CHECK(!key_absent(token));
    CHECK(entries() == std::vector<std::string>{"hdup1"});
    CHECK(segment_exists(n.segments[0]) && !segment_exists(n.segments[1]));
    check_attach_matches(first);
    CHECK(test_private_shm_destroy() == 0);
    CHECK(key_absent(token) && !segment_exists(n.segments[0]));
    std::printf("ok h10 a stage that could not create its segment removes no other\n");
}

// h11. `shm -d` that could not take the stager lock keeps the stub of a live
//      1 GB stager that has not recorded its segment (it may still create it),
//      or of one whose life cannot be checked (another PID namespace), and
//      removes it once that stager is known dead.
void check_huge_destroy_unlocked() {
    Names n({"hlive"});
    HugeCase c;
    make_stub(n.segment, 0, (uint64_t)getpid(), c.fresh(), 65536);
    register_entry("hlive", 4096);
    sem_t *lock = sem_open(private_lock, 0);   // a stager holds the lock
    CHECK(lock != SEM_FAILED);
    CHECK(sem_wait(lock) == 0);
    g_destroy_wait_ms = 100;
    const int rc = test_private_shm_destroy();
    sem_post(lock);
    sem_close(lock);
    CHECK(rc != 0);
    CHECK(segment_exists(n.segment) && entries() == std::vector<std::string>{"hlive"});
    CHECK(!lock_exists());   // the recovery unlinked the lock it could not take
    // Its PID dead here, but the stager was in another PID namespace, where
    // that PID may be a live stager: still kept.
    set_reserved(n.segment, 0, dead_pid());
    set_reserved(n.segment, 4, ipc_namespace() | (pid_namespace() + 1) << 32);
    CHECK(test_private_shm_destroy() != 0);
    CHECK(segment_exists(n.segment) && entries() == std::vector<std::string>{"hlive"});
    // Dead, seen from its own PID namespace: removed.
    set_reserved(n.segment, 4, stub_namespaces());
    CHECK(test_private_shm_destroy() == 0);
    CHECK(!segment_exists(n.segment));
    CHECK(shm_open(private_control, O_RDONLY, 0) < 0 && errno == ENOENT);
    std::printf("ok h11 destroy without the lock keeps a live 1 GB stager's stub\n");
}

// h12. A segment whose header cannot be read may be a 1 GB stub, so `shm -d`
//      keeps it and its entry rather than unlinking it blind. (Root reads
//      every segment.)
void check_unreadable_header() {
    if (geteuid() == 0) {
        std::printf("SKIP: h12 unreadable header (running as root)\n");
        return;
    }
    Names n({"hidden"});
    HugeCase c;
    const int fd = shm_open(n.segment.c_str(), O_CREAT | O_EXCL | O_RDWR, 0);
    CHECK(fd >= 0);
    CHECK(ftruncate(fd, (off_t)sizeof(bwa_shm_header_t)) == 0);
    close(fd);
    register_entry("hidden", 4096);
    CHECK(test_private_shm_destroy() != 0);
    CHECK(entries() == std::vector<std::string>{"hidden"});
    CHECK(shm_unlink(n.segment.c_str()) == 0);   // it was still there
    CHECK(test_private_shm_destroy() == 0);
    std::printf("ok h12 destroy keeps a segment whose header cannot be read\n");
}

// Run `bwa-mem3 shm <args>` (main_shm) in a child; returns its exit status and
// sets *err to what it wrote to stderr.
int run_shm(const Names &n, const std::vector<std::string> &args, std::string *err) {
    const std::string path = n.dir.file("stderr");
    fflush(NULL);
    const pid_t pid = fork();
    CHECK(pid >= 0);
    if (pid == 0) {
        alarm(20);
        if (freopen(path.c_str(), "w", stderr) == NULL) _exit(94);
        std::vector<char *> argv = {(char *)"shm"};
        for (const std::string &a : args) argv.push_back((char *)a.c_str());
        argv.push_back(NULL);
        // exit, as bwa-mem3 does: a completed 1 GB stage must survive the
        // handler that removes the segment of a stage that exits mid-pack.
        exit(main_shm((int)argv.size() - 1, argv.data()));
    }
    g_children.push_back(pid);
    const int status = wait_status(pid);
    CHECK(WIFEXITED(status));
    *err = slurp(path);
    return WEXITSTATUS(status);
}

// h9. The command line: --huge-pages reaches the stage, and is rejected with
//     -l and -d.
void check_huge_cli() {
    Names n({"hcli"});
    HugeCase c;
    const std::string prefix = n.prefix("hcli");
    std::string err;
    for (const char *flag : {"-l", "-d"}) {
        CHECK(run_shm(n, {"--huge-pages", flag}, &err) == 1);
        CHECK(has(err, "--huge-pages cannot be combined with -l or -d"));
    }
    const uint64_t token = c.fresh();
    CHECK(run_shm(n, {"--huge-pages", prefix}, &err) == 0);
    CHECK(has(err, "x 1 GB pages (System V key 0x"));
    CHECK(read_header(n.segment).magic == BWA_SHM_HUGE_MAGIC && !key_absent(token));
    CHECK(run_shm(n, {"-d"}, &err) == 0);
    CHECK(key_absent(token) && !segment_exists(n.segment));
    CHECK(run_shm(n, {prefix}, &err) == 0);
    CHECK(read_header(n.segment).magic == BWA_SHM_MAGIC);
    CHECK(run_shm(n, {"-d"}, &err) == 0);
    std::printf("ok h9 shm --huge-pages on the command line\n");
}

#ifdef __linux__
// 1 GB pages neither in use nor reserved; -1 if the host has no 1 GB pool.
long pool_available() {
    long v[2] = {0, 0};
    const char *files[2] = {"free_hugepages", "resv_hugepages"};
    for (int i = 0; i < 2; ++i) {
        const std::string path =
            std::string("/sys/kernel/mm/hugepages/hugepages-1048576kB/") + files[i];
        FILE *f = fopen(path.c_str(), "r");
        if (f == NULL) return -1;
        const int got = fscanf(f, "%ld", &v[i]);
        fclose(f);
        if (got != 1) return -1;
    }
    return v[0] - v[1];
}

// The kernel page size, in kB, of this process's mapping at `addr`.
long kernel_page_kb(const void *addr) {
    FILE *f = fopen("/proc/self/smaps", "r");
    CHECK(f != NULL);
    char line[512];
    unsigned long lo = 0, hi = 0;
    long kb = -1, found = -1;
    bool in = false;
    while (fgets(line, sizeof line, f) != NULL) {
        if (sscanf(line, "%lx-%lx ", &lo, &hi) == 2 && strchr(line, '-') < strchr(line, ' '))
            in = (unsigned long)addr >= lo && (unsigned long)addr < hi;
        else if (in && sscanf(line, "KernelPageSize: %ld kB", &kb) == 1) found = kb;
    }
    fclose(f);
    return found;
}
#endif

// hp. Real 1 GB pages (BWA_TEST_HUGETLB_POOL=1 on a Linux host with a free page in
//     the 1 GB pool and SHM_HUGETLB permission; otherwise skipped, saying
//     why): the attached index is mapped on them, and its segment is gone
//     after destroy and after the reclaim of a killed stage. The pool's counts
//     are host-wide, so they are checked only with BWA_TEST_HUGETLB_POOL=exclusive,
//     which says that nothing else on the host uses the 1 GB pool meanwhile:
//     every page is then back in the pool.
void check_huge_real_pool() {
    std::string why = "not Linux";
#ifdef __linux__
    const char *on = getenv("BWA_TEST_HUGETLB_POOL");
    const bool exclusive = on != NULL && strcmp(on, "exclusive") == 0;
    const long baseline = pool_available();
    why.clear();
    if (on == NULL || (strcmp(on, "1") != 0 && !exclusive)) {
        why = "BWA_TEST_HUGETLB_POOL=1 (or =exclusive) is not set";
    } else if (baseline < 1) {
        why = "the host's 1 GB pool has no free page";
    } else {
        const int probe = shmget(IPC_PRIVATE, (size_t)1 << 30,
                                 IPC_CREAT | SHM_HUGETLB | SHM_HUGE_1GB | 0600);
        if (probe < 0)
            why = std::string("a 1 GB System V segment cannot be created: ") + strerror(errno)
                  + " (it needs CAP_IPC_LOCK or vm.hugetlb_shm_group)";
        else
            shmctl(probe, IPC_RMID, NULL);
    }
    if (why.empty()) {
        // The pool holds `want` fewer pages than at the start (exclusive only).
        const auto pool_is = [&](long want) { CHECK(!exclusive || pool_available() == want); };
        {
            Names n({"hpool"});
            HugeCase c;
            g_huge_flags = SHM_HUGETLB | SHM_HUGE_1GB;
            g_huge_page = 1ull << 30;
            const std::string prefix = n.prefix("hpool");
            const uint64_t staged = c.fresh();
            CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, /*huge=*/true)), 0));
            pool_is(baseline - 1);
            size_t len = 0, map_len = 0;
            uint8_t *m = bwa_shm_attach(prefix.c_str(), &len, &map_len);
            CHECK(m != NULL && map_len == (1ull << 30) && kernel_page_kb(m) == 1048576);
            CHECK(munmap(m, map_len) == 0);
            check_attach_matches(prefix, (size_t)1 << 30);
            CHECK(test_private_shm_destroy() == 0);
            CHECK(key_absent(staged));
            pool_is(baseline);
            const uint64_t killed = c.fresh();
            CHECK(killed_huge_stage(wait_status(spawn_stage(prefix, HOOK_SIGKILL, true)),
                                    n.segment));
            CHECK(!key_absent(killed));
            pool_is(baseline - 1);   // held by the dead stage
            const uint64_t restaged = c.fresh();
            CHECK(exited(wait_status(spawn_stage(prefix, HOOK_NONE, true)), 0));
            CHECK(key_absent(killed) && !key_absent(restaged));
            pool_is(baseline - 1);   // reclaimed, then staged
            CHECK(test_private_shm_destroy() == 0);
            CHECK(key_absent(restaged));
        }
        pool_is(baseline);
        std::printf("ok hp real 1 GB pages: mapped on them, removed after destroy and "
                    "reclaim%s\n", exclusive ? "; every page back in the pool" : "");
        return;
    }
#endif
    std::printf("SKIP: hp real 1 GB pages (%s)\n", why.c_str());
}

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
        if (rmid_forced() != 0 && !(rmid_forced() < 0 && !is_linux())) {
            // The kernel removes every System V segment when its stager exits,
            // or may (the setting cannot be read).
            std::printf("SKIP: h1-h12, hp 1 GB staging (kernel.shm_rmid_forced %s on this "
                        "host)\n", rmid_forced() > 0 ? "=1" : "cannot be read");
        } else {
            check_huge_stage();
            check_huge_interrupted();
            check_huge_foreign();
            check_huge_rollback();
            check_huge_refused();
            check_huge_removal_fails();
            check_huge_attach_validation();
            check_huge_preflight();
            check_huge_preflight_last_page();
            check_huge_rollback_same_key();
            check_huge_destroy_unlocked();
            check_unreadable_header();
            check_huge_cli();
            check_huge_real_pool();
        }
        std::printf("shm_stage_test: OK\n");
        return 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "%s (errno=%d %s)\n", e.what(), errno, std::strerror(errno));
        return 1;
    }
}
