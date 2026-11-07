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

// Attach succeeds and its bytes equal bwa_shm_pack_from_disk's.
void check_attach_matches(const std::string &prefix) {
    size_t len = 0;
    uint8_t *m = bwa_shm_attach(prefix.c_str(), &len);
    CHECK(m != NULL);
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
    munmap(m, len);
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
pid_t spawn_stage(const std::string &prefix, HookMode mode) {
    fflush(NULL);
    const pid_t pid = fork();
    CHECK(pid >= 0);
    if (pid == 0) {
        for (int fd : g_child_close) close(fd);
        alarm(20);
        g_mode = mode;
        bwa_shm_stage_map_hook = stage_hook;
        const int rc = test_private_shm_stage(prefix.c_str());
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
    explicit PausedStager(const std::string &prefix) {
        CHECK(pipe(ready) == 0 && pipe(go) == 0);
        g_ready_fd = ready[1];
        g_go_fd = go[0];
        g_child_close.push_back(ready[0]);
        g_child_close.push_back(go[1]);
        pid = spawn_stage(prefix, HOOK_PAUSE);
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
        std::printf("shm_stage_test: OK\n");
        return 0;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "%s (errno=%d %s)\n", e.what(), errno, std::strerror(errno));
        return 1;
    }
}
