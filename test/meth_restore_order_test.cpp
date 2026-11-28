// Refuse inconsistent meth metadata before attempting a reference-sized lookup
// allocation. Returning NULL is not enough: the lookup builder treats OOM as a
// normal fallback, so the interposer exits on the attempted allocation itself.
#include "bntseq.h"
#include "bwa_shm.h"
#include "fm_index_writer.h"
#include "read_index_ele.h"
#include "framework/temp_dir.h"

#include <cerrno>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <string>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

int main_mem(int argc, char **argv);

static bool reject_large_calloc = false;
extern "C" void *calloc(size_t count, size_t size)
{
    if (size && count > SIZE_MAX / size) { errno = ENOMEM; return NULL; }
    const size_t bytes = count * size;
    if (reject_large_calloc && bytes > (1U << 20)) {
        static const char message[] = "FAIL: attempted large calloc before metadata rejection\n";
        (void)write(STDERR_FILENO, message, sizeof(message) - 1);
        _exit(90);
    }
    void *p = malloc(bytes ? bytes : 1);
    if (p) memset(p, 0, bytes);
    return p;
}

namespace {

void metadata(const std::string &prefix, int64_t length, bool seed)
{
    const int count = seed ? 2 : 1;
    std::ofstream ann(prefix + ".ann"), amb(prefix + ".amb");
    ann << length << ' ' << count << " 11\n";
    if (seed) ann << "0 rchrA (null)\n0 8 0\n0 fchrA (null)\n8 8 0\n";
    else ann << "0 chrA (null)\n0 8 0\n";
    amb << length << ' ' << count << " 0\n";
}

void small_index(const std::string &prefix, bool seed)
{
    metadata(prefix, seed ? 16 : 8, seed);
    std::ofstream pac(prefix + ".pac", std::ios::binary);
    const std::string packed(seed ? 6 : 4, '\0');
    pac.write(packed.data(), packed.size());
    const int64_t rows = seed ? 33 : 17, rate = 3;
    const Bwt2bitLayout layout = fmi_bwt2bit_layout(rows, rate);
    std::vector<char> data((size_t)layout.total_bytes, 0);
    memcpy(data.data(), &rows, sizeof rows);
    memcpy(data.data() + layout.off_sa_compx_tag, &rate, sizeof rate);
    std::ofstream bwt(prefix + ".bwt.2bit.64", std::ios::binary);
    bwt.write(data.data(), data.size());
}

enum class Load { Stage, Mem, Pac };

bool refuses(const std::string &prefix, const std::string &err, Load load,
             const char *expected)
{
    fflush(NULL);
    const pid_t pid = fork();
    if (pid == 0) {
        if (!freopen(err.c_str(), "w", stderr)) _exit(91);
        setvbuf(stderr, NULL, _IONBF, 0);
        const struct rlimit no_core = {0, 0};
        const struct rlimit data_cap = {512U << 20, 512U << 20};
        if (setrlimit(RLIMIT_CORE, &no_core)) _exit(92);
        // A backstop only: the calloc interposer is the check. Some kernels
        // (recent macOS) reject any finite RLIMIT_DATA with EINVAL.
        if (setrlimit(RLIMIT_DATA, &data_cap) && errno != EINVAL) _exit(92);
        alarm(20);
        setenv("BWA3_BWTSEED_JUMP", "0", 1);
        reject_large_calloc = true;
        if (load == Load::Stage) {
            char cmd[] = "shm", meth[] = "--meth";
            char *args[] = {cmd, meth, const_cast<char *>(prefix.c_str()), NULL};
            _exit(main_shm(3, args));
        } else if (load == Load::Mem) {
            char cmd[] = "mem", meth[] = "--meth", threads[] = "-t", one[] = "1";
            const std::string reads = prefix + ".fq";
            char *args[] = {cmd, meth, threads, one, const_cast<char *>(prefix.c_str()),
                            const_cast<char *>(reads.c_str()), NULL};
            _exit(main_mem(6, args));
        } else {
            indexEle idx;
            idx.bwa_idx_load_ele(prefix.c_str(), BWA_IDX_BNS | BWA_IDX_PAC, 1);
            _exit(0);
        }
    }
    if (pid < 0) return false;
    int status;
    if (waitpid(pid, &status, 0) != pid) return false;
    std::ifstream input(err);
    const std::string message((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 1 || message.find(expected) == std::string::npos) {
        fprintf(stderr, "FAIL: expected '%s' without large allocation (wait status %d)\n%s",
                expected, status, message.c_str());
        return false;
    }
    return true;
}

} // namespace

int main()
{
    bwa_tests::TempDir dir("meth_restore_order");
    const std::string prefix = dir.file("x" + dir.path().substr(dir.path().size() - 6));
    const std::string seed = prefix + ".meth", err = dir.file("stderr");
    std::ofstream reads(prefix + ".fq");
    reads << "@r\nACGT\n+\nIIII\n";
    reads.close();
    bool ok = true;
    for (const int64_t length : {1LL << 48, 1LL << 36}) {
        small_index(prefix, false);
        small_index(seed, true);
        metadata(prefix, length, false);
        ok = refuses(prefix, err, Load::Stage, "inconsistent reference length") && ok;
        ok = refuses(prefix, err, Load::Mem, "inconsistent reference length") && ok;
        ok = refuses(prefix, err, Load::Pac, ".pac is") && ok;
        small_index(prefix, false);
        metadata(seed, length, true);
        ok = refuses(prefix, err, Load::Stage, "inconsistent reference length") && ok;
    }
    // Matching header lengths alone are insufficient: the seed FM file is small.
    metadata(prefix, 1LL << 35, false);
    metadata(seed, 1LL << 36, true);
    ok = refuses(prefix, err, Load::Stage, "do not belong together") && ok;
    if (ok) puts("PASS: meth metadata rejects inconsistent lengths before large allocations");
    return ok ? 0 : 1;
}
