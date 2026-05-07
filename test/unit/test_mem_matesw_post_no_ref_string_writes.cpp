// Regression test: mem_matesw_batch_post must not write to mmc->ref_string.
//
// Bug: PR #76 wired mem_matesw_batch_pre/post to bns_fetch_seq_v2, which
// returns a pointer into mmc->ref_string instead of a malloc'd copy.
// mem_matesw_batch_post's `index == -1` fallback then passed that pointer
// directly to ksw_align2, whose internal revseq() reverses its `target`
// argument in place via xor-swap. With shm-backed ref_string (PROT_READ
// mmap'd /dev/shm/bwaidx-*), that in-place write SIGSEGVs.
//
// We can't detect the mutation by comparing bytes before/after the call:
// ksw_align2 calls revseq twice (lines 375 and 381 of ksw.cpp) so the net
// effect on the buffer is zero. To detect the *write attempt* itself we
// place ref_string in an mprotect(PROT_READ) page, fork, and run the
// function in the child. Pre-fix, the child SIGSEGVs at the first revseq.
// Post-fix, it exits cleanly.

#include "doctest/doctest.h"

extern "C" {
#include "bntseq.h"
}
#include "bwamem.h"

#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>

#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>

#define _set_pac(pac, l, c) ((pac)[(l) >> 2] |= (c) << ((~(l) & 3) << 1))

namespace {

// In-child callable: builds minimal bns + ms inputs, then calls
// mem_matesw_batch_post against a PROT_READ ref_string. _exit on success;
// SIGSEGV terminates if the function writes to ref_string.
[[noreturn]] void child_run_matesw_post(uint8_t *ro_ref_string,
                                        int64_t l_pac) {
    // Restore default SIGSEGV disposition so the kernel terminates the
    // child with WTERMSIG=SIGSEGV rather than running doctest's inherited
    // crash handler (which would emit a redundant doctest report).
    std::signal(SIGSEGV, SIG_DFL);
    std::signal(SIGABRT, SIG_DFL);

    // Pack the same bases that ro_ref_string holds (forward half).
    std::vector<uint8_t> pac(static_cast<size_t>((l_pac + 3) / 4 + 1), 0);
    for (int64_t i = 0; i < l_pac; ++i) {
        _set_pac(pac.data(), i, ro_ref_string[i]);
    }

    bntann1_t anns[1];
    std::memset(anns, 0, sizeof(anns));
    anns[0].offset = 0;
    anns[0].len    = static_cast<int32_t>(l_pac);
    anns[0].name   = const_cast<char *>("chr_test");
    anns[0].anno   = const_cast<char *>("");

    bntseq_t bns;
    std::memset(&bns, 0, sizeof(bns));
    bns.l_pac   = l_pac;
    bns.n_seqs  = 1;
    bns.anns    = anns;
    bns.n_holes = 0;
    bns.ambs    = nullptr;

    mem_opt_t *opt = mem_opt_init();
    if (opt == nullptr) _exit(2);

    // FF orientation only (r==0). FF skips the reverse-complement detour
    // (is_rev=false) so ms is passed straight through to ksw_align2 as
    // `query`. ksw_align2 calls revseq on `target` only when its first SW
    // score >= xtra & 0xffff (= opt->min_seed_len * opt->a = 19). Aligning
    // ms against itself (ms = first 100 bytes of the rescue region) gives
    // score 100, which clears that threshold.
    mem_pestat_t pes[4];
    std::memset(pes, 0, sizeof(pes));
    pes[1].failed = 1;
    pes[2].failed = 1;
    pes[3].failed = 1;
    pes[0].failed = 0;
    pes[0].low    = 50;
    pes[0].high   = 250;
    pes[0].avg    = 150.0;
    pes[0].std    = 30.0;

    // Anchor: forward strand at 200..300. For r=0 (is_rev=false,
    // is_larger=true): rb = a.rb + low = 250, re = a.rb + high + l_ms
    // = 550. re-rb = 300 >= opt->min_seed_len = 19.
    mem_alnreg_t a;
    std::memset(&a, 0, sizeof(a));
    a.rid           = 0;
    a.rb            = 200;
    a.re            = 300;
    a.qb            = 0;
    a.qe            = 100;
    a.score         = 100;
    a.is_alt        = 0;
    a.secondary     = -1;
    a.secondary_all = -1;

    constexpr int l_ms = 100;
    std::vector<uint8_t> ms(l_ms);
    for (int i = 0; i < l_ms; ++i) ms[i] = ro_ref_string[250 + i];

    mem_alnreg_v ma;
    std::memset(&ma, 0, sizeof(ma));

    int32_t gar[4] = {-1, -1, -1, -1};

    kswr_t *aln_storage = nullptr;
    kswr_t **myaln = &aln_storage;

    mem_cache mmc;
    std::memset(&mmc, 0, sizeof(mmc));
    mmc.ref_string = ro_ref_string;

    // The buggy index == -1 path inside this call passes a slice of
    // mmc->ref_string straight to ksw_align2, which writes via revseq.
    (void)mem_matesw_batch_post(opt, &bns, pac.data(), pes,
                                &a, l_ms, ms.data(),
                                &ma, myaln, /*gcnt*/0, gar, &mmc);

    free(opt);
    if (ma.a) free(ma.a);
    _exit(0);
}

}  // namespace

TEST_CASE("mem_matesw_batch_post does not write to read-only ref_string "
          "(shm + ksw_align2 revseq regression)") {
    constexpr int64_t L = 1024;
    const size_t page = static_cast<size_t>(sysconf(_SC_PAGESIZE));
    const size_t alloc =
        ((static_cast<size_t>(2 * L) + page - 1) / page) * page;

    // mmap a writable region, populate the ref bytes, then mprotect to
    // PROT_READ. Any write attempt by mem_matesw_batch_post downstream
    // generates SIGSEGV in the child.
    uint8_t *ref_string = static_cast<uint8_t *>(
        mmap(nullptr, alloc, PROT_READ | PROT_WRITE,
             MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
    REQUIRE(ref_string != MAP_FAILED);

    for (int64_t i = 0; i < L; ++i) {
        ref_string[i] = static_cast<uint8_t>(i & 3);
    }
    for (int64_t i = 0; i < L; ++i) {
        ref_string[2 * L - 1 - i] = static_cast<uint8_t>(3 ^ ref_string[i]);
    }
    REQUIRE(mprotect(ref_string, alloc, PROT_READ) == 0);

    pid_t pid = fork();
    REQUIRE(pid >= 0);

    if (pid == 0) {
        child_run_matesw_post(ref_string, L);
    }

    int status = 0;
    REQUIRE(waitpid(pid, &status, 0) == pid);

    munmap(ref_string, alloc);

    if (WIFSIGNALED(status)) {
        MESSAGE("child terminated by signal ", WTERMSIG(status),
                " (SIGSEGV=", SIGSEGV, ")");
    } else if (WIFEXITED(status)) {
        MESSAGE("child exited with code ", WEXITSTATUS(status));
    }

    // Pre-fix: SIGSEGV at the first revseq inside ksw_align2.
    // Post-fix: clean exit because mem_matesw_batch_post copied the slice.
    CHECK_FALSE(WIFSIGNALED(status));
    CHECK(WIFEXITED(status));
    CHECK(WEXITSTATUS(status) == 0);
}
