// Exercise real disk/private-stage loads and main_mem, without touching the
// host's shared-memory registry. Missing phases, wrong verbosity guards, and
// stale huge-page worker limits must fail independently of elapsed values.
#define bwa_shm_stage test_timing_shm_stage
#define bwa_shm_destroy test_timing_shm_destroy
#include "bwa_shm.h"
#include "framework/temp_index.h"

#include <cmath>
#include <fstream>
#include <map>
#include <semaphore.h>
#include <sstream>
#include <sys/mman.h>
#include <sys/wait.h>

static char private_control[32], private_lock[32], private_index[32];
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

int main_mem(int argc, char **argv);

namespace {

void require(bool ok, const std::string &message)
{
    if (!ok) throw std::runtime_error(message);
}

// Reserve both private names exclusively, as the existing stage regressions do.
// Destruction can only destroy this fixture's registry, never the global one.
struct PrivateStage {
    bwa_tests::TempDir dir{"lt"};
    std::string prefix;
    bool control_created = false, lock_created = false;
    PrivateStage(const std::string &source)
    {
        const std::string token = dir.path().substr(dir.path().size() - 6);
        snprintf(private_control, sizeof private_control, "/lt%sc", token.c_str());
        snprintf(private_lock, sizeof private_lock, "/lt%sl", token.c_str());
        snprintf(private_index, sizeof private_index, "/lt%s-", token.c_str());
        prefix = dir.file("x" + token); // short enough for macOS POSIX shm names
        for (const char *suffix : {".ann", ".amb", ".pac", ".bwt.2bit.64"})
            require(symlink((source + suffix).c_str(), (prefix + suffix).c_str()) == 0,
                    "cannot link timing fixture index");
        const int fd = shm_open(private_control, O_CREAT | O_EXCL | O_RDWR, 0600);
        require(fd >= 0, "cannot reserve private timing registry");
        close(fd);
        control_created = true;
        sem_t *lock = sem_open(private_lock, O_CREAT | O_EXCL, 0600, 1);
        if (lock == SEM_FAILED) {
            shm_unlink(private_control);
            throw std::runtime_error("cannot reserve private timing lock");
        }
        sem_close(lock);
        lock_created = true;
    }
    ~PrivateStage()
    {
        if (control_created && lock_created) test_timing_shm_destroy();
        if (control_created) shm_unlink(private_control);
        if (lock_created) sem_unlink(private_lock);
    }
};

// Every expected duration must occur exactly once, parse fully, and be finite
// and nonnegative. Merely matching a report prefix would miss dropped phases.
std::map<std::string, std::string> report(const std::string &log, const std::string &tag)
{
    const size_t at = log.find(tag);
    require(at != std::string::npos, "missing report: " + tag + "\n" + log);
    require(log.find(tag, at + tag.size()) == std::string::npos, "duplicate report: " + tag);
    const size_t end = log.find('\n', at);
    std::istringstream fields(log.substr(at + tag.size(), end - at - tag.size()));
    std::map<std::string, std::string> values;
    std::string field;
    while (fields >> field) {
        const size_t eq = field.find('=');
        require(eq != std::string::npos, "malformed timing field: " + field);
        require(values.emplace(field.substr(0, eq), field.substr(eq + 1)).second,
                "duplicate timing field: " + field);
    }
    return values;
}

void durations(const std::map<std::string, std::string> &fields,
               std::initializer_list<const char *> names)
{
    for (const char *name : names) {
        const auto it = fields.find(name);
        require(it != fields.end(), std::string("missing timing field: ") + name);
        char *end = NULL;
        const double value = strtod(it->second.c_str(), &end);
        require(end != it->second.c_str() && *end == '\0' && std::isfinite(value) && value >= 0,
                std::string("invalid duration: ") + name + "=" + it->second);
    }
}

const char *const phase_tag = "[M::load_index] phases (s):";
const char *const wall_tag = "[M::main_mem] wall (s):";

std::string run_mem(const PrivateStage &stage, int verbosity, const std::string &stem,
                    bool jump = false)
{
    const std::string sam = stage.dir.file(stem + ".sam"), err = stage.dir.file(stem + ".err");
    const std::string reads = stage.dir.file("reads.fq");
    fflush(NULL);
    const pid_t pid = fork();
    require(pid >= 0, "cannot fork mem timing check");
    if (pid == 0) {
        alarm(20);
        if (!freopen(sam.c_str(), "w", stdout) || !freopen(err.c_str(), "w", stderr)) _exit(90);
        setenv("BWA3_BWTSEED_JUMP", jump ? "4" : "0", 1);
        unsetenv("BWA3_LOAD_THREADS");
        char cmd[] = "mem", threads[] = "-t", two[] = "2", verbose[] = "-v", four[] = "4";
        char *args[] = {cmd, threads, two, verbose, four,
                        const_cast<char *>(stage.prefix.c_str()), const_cast<char *>(reads.c_str()), NULL};
        char *defaults[] = {cmd, threads, two, const_cast<char *>(stage.prefix.c_str()),
                            const_cast<char *>(reads.c_str()), NULL};
        const int rc = verbosity == 4 ? main_mem(7, args) : main_mem(5, defaults);
        fflush(NULL);
        _exit(rc);
    }
    int status = 0;
    require(waitpid(pid, &status, 0) == pid, "cannot wait for mem timing check");
    const std::string log = stage.dir.read(stem + ".err");
    require(WIFEXITED(status) && WEXITSTATUS(status) == 0,
            "mem timing check failed (wait status " + std::to_string(status) + ")\n" + log);
    if (verbosity == 4) {
        const auto phases = report(log, phase_tag);
        const bool staged = stem.find("shm") == 0;
        require(phases.at("source") == (staged ? "shm" : "disk"), "wrong load timing source");
        require(phases.at("jump_threads") == (jump ? "2" : "0"), "wrong jump-table worker limit");
        if (staged)
            durations(phases, {"shm_attach", "shm_fmi", "shm_bns_pac", "jump_table", "total"});
        else
            durations(phases, {"shm_attach", "open_header", "cp_occ", "sa_hi_read_pack",
                               "sa_ls", "jump_table", "bns_pac", "total"});
        durations(report(log, wall_tag), {"setup", "index", "pre_pipeline", "pipeline", "teardown"});
    } else {
        require(log.find(phase_tag) == std::string::npos && log.find(wall_tag) == std::string::npos,
                "default verbosity printed phase/wall timings\n" + log);
    }
    std::istringstream input(stage.dir.read(stem + ".sam"));
    std::string line, records;
    while (std::getline(input, line))
        if (!line.empty() && line[0] != '@') records += line + '\n';
    require(!records.empty(), "mem timing check emitted no alignment records");
    return records;
}

void check_worker_limits(const PrivateStage &stage)
{
    // The existing plain-mapping seam runs successful allocation without a
    // host huge-page pool. Fallback is unconditional and portable.
    for (bool compact : {false, true}) {
        for (const char *override_value : {static_cast<const char *>(NULL), "", "abc", "3"}) {
            const bool overridden = override_value != NULL && strcmp(override_value, "3") == 0;
            for (bwamem_hugetlb_mode_t mode : {BWAMEM_HUGETLB_FORCE_FAIL, BWAMEM_HUGETLB_TEST_PLAIN}) {
                if (override_value) setenv("BWA3_LOAD_THREADS", override_value, 1);
                else unsetenv("BWA3_LOAD_THREADS");
                const int ordinary = overridden ? 3 : index_load_threads_host_cap();
                int occurrence = ordinary;
#if defined(__linux__)
                if (!overridden && mode == BWAMEM_HUGETLB_TEST_PLAIN) occurrence = 16;
#endif
                std::string log;
                {
                    bwa_tests::StderrCapture capture;
                    FMI_search fmi(stage.prefix.c_str());
                    fmi.set_compact_occ(compact);
                    fmi.load_index(true, 16, false, mode);
                    log = capture.text();
                }
                const std::string wanted = "index load worker limits: ordinary=" + std::to_string(ordinary)
                    + " occurrence=" + std::to_string(occurrence) + " jump=0";
                require(log.find(wanted) != std::string::npos, "wrong worker diagnostic\n" + log);
                const auto phases = report(log, phase_tag);
                require(phases.at("ordinary_threads") == std::to_string(ordinary), "wrong ordinary limit");
                require(phases.at("occurrence_threads") == std::to_string(occurrence), "wrong occurrence limit");
                require(phases.at("jump_threads") == "0", "disabled jump table reported workers");
            }
        }
    }
    unsetenv("BWA3_LOAD_THREADS");
}

} // namespace

int main()
{
    try {
        unsetenv("BWA3_LOAD_THREADS");
        unsetenv("BWA3_COMPACT_OCC");
        setenv("BWA3_BWTSEED_JUMP", "0", 1);
        setenv("OMP_NUM_THREADS", "2", 1);
        std::mt19937 rng(20261125);
        const std::string sequence = bwa_tests::random_acgt(rng, 2200);
        bwa_tests::TempIndex index({{"chrA", sequence}});
        PrivateStage stage(index.prefix());
        const std::string read = sequence.substr(400, 150);
        stage.dir.write("reads.fq", "@r\n" + read + "\n+\n" + std::string(read.size(), 'I') + "\n");
        bwa_verbose = 4;
        check_worker_limits(stage);
        bwa_verbose = 3;
        const std::string disk = run_mem(stage, 3, "disk-default");
        require(run_mem(stage, 4, "disk-v4") == disk, "disk verbosity changed alignment records");
        require(run_mem(stage, 4, "disk-jump", true) == disk, "disk jump table changed alignment records");
        require(test_timing_shm_stage(stage.prefix.c_str()) == 0, "cannot stage private timing index");
        require(run_mem(stage, 3, "shm-default") == disk, "staged load changed alignment records");
        require(run_mem(stage, 4, "shm-v4") == disk, "staged verbosity changed alignment records");
        require(run_mem(stage, 4, "shm-jump", true) == disk, "staged jump table changed alignment records");
        require(test_timing_shm_destroy() == 0, "cannot remove private timing stage");
        puts("PASS: disk/private-shm phase and run-wall timings, verbosity, records and worker limits");
        return 0;
    } catch (const std::exception &e) {
        fprintf(stderr, "FAIL: index_load_timing_test: %s\n", e.what());
        return 1;
    }
}
