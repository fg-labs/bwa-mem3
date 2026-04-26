#include "libsais_build.h"

#include "fm_index_writer.h"
#include "index_prelude.h"
#include "io_utils.h"
#include "macro.h"
#include "packed_text.h"
#include "utils.h"

#include <atomic>
#include <cerrno>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

// libsais headers use `SA` as a parameter name in their public API; bwa's
// macro.h defines SA as a preprocessor constant. Undefine around the include
// so we get the real declarations. We don't need the SA macro here anyway.
#ifdef SA
# pragma push_macro("SA")
# undef SA
# define SA_MACRO_PUSHED 1
#endif
extern "C" {
#include "libsais.h"
#include "libsais64.h"
}
#ifdef SA_MACRO_PUSHED
# pragma pop_macro("SA")
# undef SA_MACRO_PUSHED
#endif

namespace {

using clock_t_ = std::chrono::steady_clock;
double elapsed_since(clock_t_::time_point t) {
    return std::chrono::duration<double>(clock_t_::now() - t).count();
}


std::string fmt_bytes(int64_t b) {
    char out[32];
    if      (b >= (1LL << 30)) std::snprintf(out, sizeof(out), "%.1f GiB", (double)b / (double)(1LL << 30));
    else if (b >= (1LL << 20)) std::snprintf(out, sizeof(out), "%.1f MiB", (double)b / (double)(1LL << 20));
    else                        std::snprintf(out, sizeof(out), "%lld B",   (long long)b);
    return std::string(out);
}

// Pack the forward-only .pac (l_pac bases) plus its reverse-complement into
// a doubled .pac on disk. The sais-lite path built this doubled text in
// memory via pac2nt; writing it packed lets PackedText / emit_0123 /
// compute_counts consume both strands without another ASCII round-trip.
// Parallelises by striping output bytes: each thread fills a disjoint
// [byte_start, byte_end) range of the output buffer, so there's no
// aliasing at the forward/reverse-complement seam even when l_pac isn't
// a multiple of 4.
void write_doubled_pac(const char* fwd_pac_path, int64_t l_pac,
                       const char* out_path, int num_threads)
{
    PackedText fwd(fwd_pac_path, l_pac);
    int64_t doubled = 2 * l_pac;
    int64_t nbytes  = (doubled + 3) / 4;
    std::vector<uint8_t> buf((size_t)nbytes, 0);

    auto base_at = [&](int64_t idx) -> uint8_t {
        // Forward half: fwd[idx]. RC half: complement of fwd[l_pac-1-(idx-l_pac)].
        if (idx < l_pac) return fwd.get_base(idx);
        int64_t src = l_pac - 1 - (idx - l_pac);
        return (uint8_t)(3u - fwd.get_base(src));
    };

#ifdef _OPENMP
    #pragma omp parallel for num_threads(num_threads) schedule(static)
#endif
    for (int64_t by = 0; by < nbytes; ++by) {
        uint8_t b = 0;
        for (int k = 0; k < 4; ++k) {
            int64_t idx = by * 4 + k;
            if (idx >= doubled) break;
            uint8_t base = base_at(idx);
            b |= (uint8_t)(base << ((3 - k) << 1));
        }
        buf[(size_t)by] = b;
    }

    int fd = open(out_path, O_WRONLY | O_CREAT | O_TRUNC, 0600);
    if (fd < 0) err_fatal(__func__, "open('%s'): %s", out_path, strerror(errno));
    // Linux caps a single write(2) on regular files at 0x7ffff000 bytes (~2 GiB)
    // regardless of the requested length; for references where doubled-pac
    // grows past that (large custom genomes / methylation doubled wheat) the
    // raw write() returned a short count and aborted with a half-written file
    // in $TMPDIR. pwrite_all loops on EINTR + short writes.
    pwrite_all(fd, buf.data(), (size_t)nbytes, 0, "doubled.pac data");
    off_t cur = (off_t)nbytes;
    uint8_t ct = 0;
    if (doubled % 4 == 0) {
        pwrite_all(fd, &ct, 1, cur, "doubled.pac pad-byte");
        ++cur;
    }
    ct = (uint8_t)(doubled % 4);
    pwrite_all(fd, &ct, 1, cur, "doubled.pac count-byte");
    if (close(fd) != 0)
        err_fatal(__func__, "close('%s'): %s", out_path, strerror(errno));
}

std::string resolve_tmpdir(const LibsaisBuildOpts& opts) {
    if (!opts.tmpdir.empty()) return opts.tmpdir;
    if (const char* e = std::getenv("TMPDIR")) return std::string(e);
    return std::string("/tmp");
}

} // anonymous namespace

int libsais_build_fm_index(const char* prefix, int64_t pac_len,
                           const LibsaisBuildOpts& opts)
{
    int     T_user = opts.num_threads > 0 ? opts.num_threads : 1;
    int64_t N = pac_len;   // doubled-text length
    int64_t l_pac = N / 2;

    // Cap effective parallelism at (N+1) / 16 Mbp. On tiny inputs the OMP
    // pool spin-up + per-thread heap arenas cost more than the parallel
    // work saves, which breaks the --max-memory contract on small fixtures.
    constexpr int64_t kMinBasesPerThread = 16LL << 20;
    int T = (int)std::min<int64_t>(T_user,
                                   std::max<int64_t>(1, (N + 1) / kMinBasesPerThread));

    auto t0 = clock_t_::now();
    std::fprintf(stderr, "[libsais_build] N=%lld threads=%d (user-requested=%d)\n",
            (long long)N, T, T_user);

    // Per-base preflight headroom: 6 B/base (int32) / 12 B/base (int64)
    // covers measured peaks on chr22 and chr1-6 with ~20% margin for
    // libsais aux arrays and OMP/mimalloc overhead on small inputs.
    const int64_t kSais32Threshold = (int64_t)INT32_MAX - 10000;
    const bool    use_int64_sa     = (N + 1) >= kSais32Threshold;
    const int64_t per_base_bytes   = use_int64_sa ? 12 : 6;
    const int64_t est_bytes        = (N + 1) * per_base_bytes;
    if (opts.max_memory_bytes > 0 && est_bytes > opts.max_memory_bytes) {
        std::fprintf(stderr,
                "ERROR: libsais build would use ~%s (%d-bit SA); --max-memory is %s.\n",
                fmt_bytes(est_bytes).c_str(),
                use_int64_sa ? 64 : 32,
                fmt_bytes(opts.max_memory_bytes).c_str());
        std::fprintf(stderr,
                "       Raise --max-memory or use a smaller reference. "
                "(Bounded-memory SA construction is not yet implemented.)\n");
        return 3;
    }
    std::fprintf(stderr, "[libsais_build] estimate ~%s (budget %s)\n",
            fmt_bytes(est_bytes).c_str(),
            opts.max_memory_bytes > 0 ? fmt_bytes(opts.max_memory_bytes).c_str() : "unset");

    // Phase 0: materialise a doubled .pac so downstream helpers can mmap it
    // via PackedText. bns_fasta2bntseq writes forward-only; we append the
    // RC half in a scratch temp file.
    std::string tmp_base = resolve_tmpdir(opts);
    char workdir[PATH_MAX];
    std::snprintf(workdir, sizeof(workdir), "%s/bwa_libsais_XXXXXX", tmp_base.c_str());
    if (!mkdtemp(workdir))
        err_fatal(__func__, "mkdtemp under '%s' failed: %s",
                  tmp_base.c_str(), strerror(errno));
    std::string doubled_pac = std::string(workdir) + "/doubled.pac";
    // RAII guard so the workdir + doubled.pac are removed on every exit
    // path, including the err_fatal calls below. err_fatal calls exit(),
    // so we also wire an atexit hook for non-stack-unwinding paths.
    struct TmpGuard {
        std::string dir;
        std::string pac;
        static TmpGuard*& current() { static TmpGuard* g = nullptr; return g; }
        static void atexit_cleanup() {
            if (auto* g = current()) g->cleanup();
        }
        void cleanup() {
            if (!pac.empty()) { ::unlink(pac.c_str()); pac.clear(); }
            if (!dir.empty()) { ::rmdir(dir.c_str());  dir.clear(); }
        }
        ~TmpGuard() { cleanup(); current() = nullptr; }
    } tmp_guard{workdir, doubled_pac};
    TmpGuard::current() = &tmp_guard;
    std::atexit(&TmpGuard::atexit_cleanup);

    write_doubled_pac((std::string(prefix) + ".pac").c_str(), l_pac,
                      doubled_pac.c_str(), T);

    PackedText pac(doubled_pac.c_str(), N);

    if (emit_0123(pac, prefix, T) != 0)
        err_fatal(__func__, "emit_0123 failed");

    int64_t count[5] = {0};
    compute_counts(pac, count, T);
    std::fprintf(stderr, "[libsais_build] phase 0 (doubled-pac + .0123 + count) %.2fs\n",
            elapsed_since(t0));

    // Phase 1: unpack .pac into a libsais-ready byte buffer. Alphabet is
    // {0=$, 1..4=ACGT}; the trailing 0 at index N is the GSA terminator.
    auto t1 = clock_t_::now();
    std::vector<uint8_t> buf((size_t)(N + 1));
    // Capture out-of-range bases without crashing inside the OpenMP region:
    // err_fatal terminates via exit(), which interacts poorly with libomp's
    // teardown when other threads are mid-loop. Record the first offending
    // index atomically and report after the parallel region.
    std::atomic<int64_t> bad_pos{-1};
    std::atomic<unsigned> bad_val{0};
#ifdef _OPENMP
    #pragma omp parallel for num_threads(T) schedule(static)
#endif
    for (int64_t i = 0; i < N; ++i) {
        uint8_t b = pac.get_base(i);
        if (b >= 4) {
            int64_t expect = -1;
            if (bad_pos.compare_exchange_strong(expect, i)) bad_val.store(b);
        }
        buf[(size_t)i] = (uint8_t)(b + 1);
    }
    buf[(size_t)N] = 0;
    if (bad_pos.load() >= 0)
        err_fatal(__func__,
                  "unexpected non-2bit base %u at doubled-pac[%lld]",
                  bad_val.load(), (long long)bad_pos.load());
    std::fprintf(stderr, "[libsais_build] phase 1 (unpack %lld bases) %.2fs\n",
            (long long)(N + 1), elapsed_since(t1));

    // Phase 2: libsais GSA. Extra-space `fs` lets libsais avoid a realloc
    // inside the induced-sorting passes; 10000 matches upstream's example.
    auto t2 = clock_t_::now();
    const int64_t fs = 10000;
    std::vector<int64_t> sa64;
    std::vector<int32_t> sa32;
    if (!use_int64_sa) {
        sa32.resize((size_t)(N + 1 + fs));
        int32_t rc = libsais_gsa_omp(buf.data(), sa32.data(),
                                     (int32_t)(N + 1), (int32_t)fs,
                                     /*freq=*/nullptr, T);
        if (rc != 0) err_fatal(__func__, "libsais_gsa_omp failed (rc=%d)", rc);
    } else {
        sa64.resize((size_t)(N + 1 + fs));
        int64_t rc = libsais64_gsa_omp(buf.data(), sa64.data(),
                                       N + 1, fs,
                                       /*freq=*/nullptr, T);
        if (rc != 0) err_fatal(__func__, "libsais64_gsa_omp failed (rc=%lld)", (long long)rc);
    }
    std::fprintf(stderr, "[libsais_build] phase 2 (libsais_gsa) %.2fs\n", elapsed_since(t2));

    // Phase 3: streaming FM-index writer reads BWT bytes and SA samples
    // in flight from sa + buf; no bwt[] or samples[] is materialised.
    auto t3 = clock_t_::now();
    std::string bwt_path = std::string(prefix) + ".bwt.2bit.64";
    int64_t sentinel_index = -1;
    // Drive the dispatch off the explicit width flag rather than
    // sa32.empty(); the latter is correct today but is load-bearing on
    // a side effect of the resize above and would silently flip if a
    // future refactor reserves both vectors for RAII reasons.
    const void* sa_ptr = use_int64_sa
        ? (const void*)sa64.data()
        : (const void*)sa32.data();
    write_fm_index_streaming(
        bwt_path.c_str(), buf.data(),
        sa_ptr,
        /*sa_is_64bit=*/use_int64_sa,
        N, count, &sentinel_index, T);
    std::fprintf(stderr,
            "[libsais_build] phase 3 (streaming FM-index write) %.2fs, sentinel_index=%lld; total %.2fs\n",
            elapsed_since(t3), (long long)sentinel_index, elapsed_since(t0));

    return 0;
}
