// Unit tests for the THP-fallback stderr check (src/bwa_thp.h): parsing of
// /proc/self/smaps, /proc/self/status and the transparent_hugepage sysfs knobs,
// the eligibility-based decision, the per-batch state machine (driven through
// an injected reader), the sampling schedule, and the message text. All inputs
// are synthetic; the real /proc reads are a thin wrapper validated on Linux.

#include "doctest/doctest.h"

#include "bwa_thp.h"

#include <algorithm>
#include <cstring>
#include <map>
#include <string>
#include <vector>

static const long long GiB_kB = 1024LL * 1024LL;

// One VMA record in the shape of a 6.1 kernel /proc/self/smaps entry.
static std::string vma(const char *range, long long anon_kb, long long huge_kb, int eligible)
{
    std::string s = std::string(range) + " rw-p 00000000 00:00 0 \n";
    s += "Size:           " + std::to_string(anon_kb + 4) + " kB\n"
         "KernelPageSize:        4 kB\n"
         "MMUPageSize:           4 kB\n"
         "Rss:            " + std::to_string(anon_kb) + " kB\n"
         "Pss:            " + std::to_string(anon_kb) + " kB\n"
         "Pss_Anon:       " + std::to_string(anon_kb + 7) + " kB\n"
         "Shared_Clean:          0 kB\n"
         "Private_Dirty:  " + std::to_string(anon_kb) + " kB\n"
         "Referenced:     " + std::to_string(anon_kb) + " kB\n"
         "Anonymous:      " + std::to_string(anon_kb) + " kB\n"
         "LazyFree:              0 kB\n"
         "AnonHugePages:  " + std::to_string(huge_kb) + " kB\n"
         "ShmemPmdMapped:        0 kB\n"
         "FilePmdMapped:         0 kB\n"
         "Shared_Hugetlb:        0 kB\n"
         "Private_Hugetlb:       0 kB\n"
         "Swap:                  0 kB\n"
         "SwapPss:               0 kB\n"
         "Locked:                0 kB\n";
    if (eligible >= 0) s += "THPeligible:    " + std::to_string(eligible) + "\n";
    s += "VmFlags: rd wr mr mw me ac hg \n";
    return s;
}

static bwa_thp_usage_t parse_smaps(const std::string &text, size_t chunk = 0)
{
    bwa_thp_smaps_parser_t p;
    bwa_thp_smaps_init(&p);
    if (chunk == 0) chunk = text.size() ? text.size() : 1;
    for (size_t i = 0; i < text.size(); i += chunk)
        bwa_thp_smaps_feed(&p, text.data() + i, std::min(chunk, text.size() - i));
    bwa_thp_usage_t u;
    bwa_thp_smaps_finish(&p, &u);
    return u;
}

static bwa_thp_usage_t usage(long long anon_kb, long long huge_kb,
                             long long elig_kb = -1, long long elig_huge_kb = -1)
{
    bwa_thp_usage_t u;
    std::memset(&u, 0, sizeof u);
    u.anon_kb = anon_kb;
    u.anon_huge_kb = huge_kb;
    u.elig_kb = elig_kb < 0 ? anon_kb : elig_kb;
    u.elig_huge_kb = elig_huge_kb < 0 ? huge_kb : elig_huge_kb;
    u.vmas = 3;
    u.elig_fields = 3;
    return u;
}

static std::string status(long long rss_anon_kb, int thp_enabled)
{
    std::string s = "Name:\tbwa-mem3\nUmask:\t0022\nGroups:\t10 4 190\nVmRSS:\t29400000 kB\n";
    if (rss_anon_kb >= 0) s += "RssAnon:\t" + std::to_string(rss_anon_kb) + " kB\n";
    s += "RssFile:\t   10000 kB\nRssShmem:\t       0 kB\nThreads:\t19\n";
    if (thp_enabled >= 0) s += "THP_enabled:\t" + std::to_string(thp_enabled) + "\n";
    return s;
}

TEST_CASE("thp: smaps parsing sums Anonymous/AnonHugePages by THPeligible"
          * doctest::test_suite("unit/bwa_thp")) {
    const std::string text =
        vma("7f0000000000-7f0400000000", 10 * GiB_kB, 10 * GiB_kB, 1)       // eligible + huge
        + vma("7f0400000000-7f0600000000", 3 * GiB_kB, 1 * GiB_kB, 1)       // eligible, 2 GiB not huge
        + vma("7f0600000000-7f0800000000", 5 * GiB_kB, 0, 0)                // not eligible (unadvised)
        + vma("ffffd0000000-ffffd0021000", 132, 0, 0);                      // stack
    SUBCASE("whole text at once") {
        const bwa_thp_usage_t u = parse_smaps(text);
        CHECK(u.vmas == 4);
        CHECK(u.elig_fields == 4);
        CHECK(u.anon_kb == 18 * GiB_kB + 132);
        CHECK(u.anon_huge_kb == 11 * GiB_kB);
        CHECK(u.elig_kb == 13 * GiB_kB);
        CHECK(u.elig_huge_kb == 11 * GiB_kB);
    }
    SUBCASE("lines split across reads of any size give the same totals") {
        const bwa_thp_usage_t a = parse_smaps(text);
        for (size_t chunk : {1u, 3u, 7u, 63u, 64u, 65u, 4096u}) {
            const bwa_thp_usage_t b = parse_smaps(text, chunk);
            CHECK(b.anon_kb == a.anon_kb);
            CHECK(b.anon_huge_kb == a.anon_huge_kb);
            CHECK(b.elig_kb == a.elig_kb);
            CHECK(b.elig_huge_kb == a.elig_huge_kb);
            CHECK(b.vmas == a.vmas);
        }
    }
    SUBCASE("a header longer than the line buffer, and no trailing newline") {
        std::string t = vma("7f0000000000-7f0400000000", 6 * GiB_kB, 2 * GiB_kB, 1);
        const std::string path(500, 'x');
        t = "7f0400000000-7f0400200000 r--p 00000000 103:01 42 /" + path + "\nAnonymous: 8 kB\n" + t;
        t += "7f0500000000-7f0500200000 rw-p 00000000 00:00 0\nAnonymous: 16 kB";   // no '\n'
        const bwa_thp_usage_t u = parse_smaps(t, 50);
        CHECK(u.vmas == 3);
        CHECK(u.anon_kb == 6 * GiB_kB + 24);
        CHECK(u.elig_kb == 6 * GiB_kB);
        CHECK(u.elig_fields == 1);
    }
    SUBCASE("a kernel without THPeligible (< 5.0) reports no eligibility fields") {
        const bwa_thp_usage_t u = parse_smaps(vma("7f0000000000-7f0400000000", 8 * GiB_kB, 0, -1));
        CHECK(u.vmas == 1);
        CHECK(u.anon_kb == 8 * GiB_kB);
        CHECK(u.elig_fields == 0);
        CHECK(u.elig_kb == 0);
    }
    SUBCASE("fields are matched at line start, not inside a longer key") {
        const std::string t = "7f0000000000-7f0000200000 rw-p 00000000 00:00 0\n"
                              "XAnonymous: 5 kB\nPss_AnonHugePages: 7 kB\nAnonHugePages: 2 kB\n"
                              "Anonymous: 9 kB\nXTHPeligible: 1\nTHPeligible: 0\n";
        const bwa_thp_usage_t u = parse_smaps(t);
        CHECK(u.anon_kb == 9);
        CHECK(u.anon_huge_kb == 2);
        CHECK(u.elig_kb == 0);
        CHECK(u.elig_fields == 1);
    }
    SUBCASE("empty input") {
        const bwa_thp_usage_t u = parse_smaps("");
        CHECK(u.vmas == 0);
        CHECK(u.anon_kb == 0);
    }
}

TEST_CASE("thp: status parsing reads RssAnon and THP_enabled in full"
          * doctest::test_suite("unit/bwa_thp")) {
    bwa_thp_status_t s = bwa_thp_parse_status(status(25165824, 1).c_str());
    CHECK(s.rss_anon_kb == 25165824);
    CHECK(s.thp_enabled == 1);
    s = bwa_thp_parse_status(status(-1, 0).c_str());
    CHECK(s.rss_anon_kb == -1);
    CHECK(s.thp_enabled == 0);
    s = bwa_thp_parse_status(status(100, -1).c_str());   // kernel < 5.0
    CHECK(s.thp_enabled == -1);
    s = bwa_thp_parse_status(NULL);
    CHECK(s.rss_anon_kb == -1);
    CHECK(s.thp_enabled == -1);
    SUBCASE("a long Groups: line does not hide RssAnon") {
        std::string g = "Name:\tbwa-mem3\nGroups:\t";
        for (int i = 0; i < 2500; ++i) g += std::to_string(1000000 + i) + " ";
        g += "\n" + status(4242, 1);
        REQUIRE(g.size() > 16384);
        s = bwa_thp_parse_status(g.c_str());
        CHECK(s.rss_anon_kb == 4242);
        CHECK(s.thp_enabled == 1);
    }
    SUBCASE("a longer key containing the needle on an earlier line is skipped") {
        s = bwa_thp_parse_status("XRssAnon:\t5 kB\nRssAnon:\t9 kB\nTHP_enabled:\t7\n");
        CHECK(s.rss_anon_kb == 9);
        CHECK(s.thp_enabled == -1);   // not 0/1
    }
}

TEST_CASE("thp: incomplete walks are detected via RssAnon" * doctest::test_suite("unit/bwa_thp")) {
    CHECK_FALSE(bwa_thp_sample_complete(10 * GiB_kB, 24 * GiB_kB));
    CHECK(bwa_thp_sample_complete(24 * GiB_kB, 24 * GiB_kB));
    CHECK(bwa_thp_sample_complete(9 * GiB_kB, 10 * GiB_kB));
    CHECK_FALSE(bwa_thp_sample_complete(9 * GiB_kB - 1, 10 * GiB_kB));
    CHECK(bwa_thp_sample_complete(25 * GiB_kB, 24 * GiB_kB));   // grew between reads
    CHECK(bwa_thp_sample_complete(1, -1));                     // unknown RssAnon accepts
}

TEST_CASE("thp: enabled-mode parsing and config resolution" * doctest::test_suite("unit/bwa_thp")) {
    CHECK(bwa_thp_parse_enabled("always [madvise] never\n") == BWA_THP_MODE_MADVISE);
    CHECK(bwa_thp_parse_enabled("[always] madvise never\n") == BWA_THP_MODE_ALWAYS);
    CHECK(bwa_thp_parse_enabled("always madvise [never]\n") == BWA_THP_MODE_NEVER);
    CHECK(bwa_thp_parse_enabled("always [inherit] madvise never\n") == BWA_THP_MODE_INHERIT);
    CHECK(bwa_thp_parse_enabled("always madvise never\n") == BWA_THP_MODE_UNKNOWN);
    CHECK(bwa_thp_parse_enabled("[neverland]\n") == BWA_THP_MODE_UNKNOWN);
    CHECK(bwa_thp_parse_enabled("") == BWA_THP_MODE_UNKNOWN);
    CHECK(bwa_thp_parse_enabled(NULL) == BWA_THP_MODE_UNKNOWN);

    const bwa_thp_mode_t off[2] = { BWA_THP_MODE_NEVER, BWA_THP_MODE_NEVER };
    SUBCASE("no per-size knobs (< 6.8): the global mode") {
        bwa_thp_config_t c = bwa_thp_resolve_config(BWA_THP_MODE_MADVISE, BWA_THP_MODE_INHERIT, NULL, 0);
        CHECK(c.pmd_mode == BWA_THP_MODE_MADVISE);
        CHECK(std::string(c.knob) == BWA_THP_KNOB_GLOBAL);
        c = bwa_thp_resolve_config(BWA_THP_MODE_UNKNOWN, BWA_THP_MODE_INHERIT, NULL, 0);
        CHECK(c.pmd_mode == BWA_THP_MODE_UNKNOWN);
        CHECK(c.knob == NULL);
        c = bwa_thp_resolve_config(BWA_THP_MODE_INHERIT, BWA_THP_MODE_INHERIT, NULL, 0);
        CHECK(c.pmd_mode == BWA_THP_MODE_UNKNOWN);
    }
    SUBCASE("hugepages-2048kB overrides the global mode unless it inherits") {
        bwa_thp_config_t c = bwa_thp_resolve_config(BWA_THP_MODE_MADVISE, BWA_THP_MODE_NEVER, off, 2);
        CHECK(c.pmd_mode == BWA_THP_MODE_NEVER);
        CHECK(std::string(c.knob) == BWA_THP_KNOB_PMD);
        c = bwa_thp_resolve_config(BWA_THP_MODE_NEVER, BWA_THP_MODE_ALWAYS, off, 2);
        CHECK(c.pmd_mode == BWA_THP_MODE_ALWAYS);
        CHECK(std::string(c.knob) == BWA_THP_KNOB_PMD);
        c = bwa_thp_resolve_config(BWA_THP_MODE_NEVER, BWA_THP_MODE_INHERIT, off, 2);
        CHECK(c.pmd_mode == BWA_THP_MODE_NEVER);
        CHECK(std::string(c.knob) == BWA_THP_KNOB_GLOBAL);
        c = bwa_thp_resolve_config(BWA_THP_MODE_ALWAYS, BWA_THP_MODE_UNKNOWN, off, 2);
        CHECK(c.pmd_mode == BWA_THP_MODE_UNKNOWN);
    }
    SUBCASE("another mTHP size enabled puts the run out of scope") {
        const bwa_thp_mode_t on[2]  = { BWA_THP_MODE_NEVER, BWA_THP_MODE_ALWAYS };
        const bwa_thp_mode_t inh[2] = { BWA_THP_MODE_INHERIT, BWA_THP_MODE_NEVER };
        const bwa_thp_mode_t unk[2] = { BWA_THP_MODE_NEVER, BWA_THP_MODE_UNKNOWN };
        // an unparseable per-size knob is out of scope too
        CHECK(bwa_thp_resolve_config(BWA_THP_MODE_MADVISE, BWA_THP_MODE_INHERIT, unk, 2).pmd_mode
              == BWA_THP_MODE_UNKNOWN);
        CHECK(bwa_thp_resolve_config(BWA_THP_MODE_MADVISE, BWA_THP_MODE_INHERIT, on, 2).pmd_mode
              == BWA_THP_MODE_UNKNOWN);
        CHECK(bwa_thp_resolve_config(BWA_THP_MODE_MADVISE, BWA_THP_MODE_INHERIT, inh, 2).pmd_mode
              == BWA_THP_MODE_UNKNOWN);
        // inherit of a global `never` is off, so 2 MB alone is judged.
        CHECK(bwa_thp_resolve_config(BWA_THP_MODE_NEVER, BWA_THP_MODE_ALWAYS, inh, 2).pmd_mode
              == BWA_THP_MODE_ALWAYS);
    }
}

TEST_CASE("thp: decision thresholds" * doctest::test_suite("unit/bwa_thp")) {
    const bwa_thp_mode_t MA = BWA_THP_MODE_MADVISE, AL = BWA_THP_MODE_ALWAYS, NV = BWA_THP_MODE_NEVER;
    SUBCASE("fully-THP human WGS runs never warn (<= ~0.3 GiB short of ~27.6)") {
        bwa_thp_usage_t a = usage(27600 * 1024LL, 27300 * 1024LL);
        bwa_thp_usage_t b = usage(28 * GiB_kB, 28 * GiB_kB);
        CHECK(bwa_thp_assess(MA, 1, &a) == BWA_THP_OK);
        CHECK(bwa_thp_assess(AL, 1, &b) == BWA_THP_OK);
    }
    SUBCASE("fragmented runs warn (1.4-5 GiB short of ~27.6)") {
        bwa_thp_usage_t a = usage(27600 * 1024LL, 26200 * 1024LL);
        bwa_thp_usage_t b = usage(27600 * 1024LL, 22600 * 1024LL);
        CHECK(bwa_thp_assess(MA, 1, &a) == BWA_THP_WARN_FALLBACK);
        CHECK(bwa_thp_assess(AL, -1, &b) == BWA_THP_WARN_FALLBACK);
    }
    SUBCASE("only eligible memory counts: unadvised 4 KB memory is not a fallback") {
        // 28 GiB anonymous, 10 GiB of it eligible and fully huge, 18 GiB unadvised.
        bwa_thp_usage_t a = usage(28 * GiB_kB, 10 * GiB_kB, 10 * GiB_kB, 10 * GiB_kB);
        CHECK(bwa_thp_assess(MA, 1, &a) == BWA_THP_OK);
        // the same with 2 GiB of the eligible 10 not huge: a fallback.
        bwa_thp_usage_t b = usage(28 * GiB_kB, 8 * GiB_kB, 10 * GiB_kB, 8 * GiB_kB);
        CHECK(bwa_thp_assess(MA, 1, &b) == BWA_THP_WARN_FALLBACK);
        // no THPeligible lines at all (< 5.0): never a fallback verdict.
        bwa_thp_usage_t c = usage(28 * GiB_kB, 0, 0, 0);
        c.elig_fields = 0;
        CHECK(bwa_thp_assess(MA, 1, &c) == BWA_THP_OK);
        CHECK(bwa_thp_assess(NV, 1, &c) == BWA_THP_WARN_DISABLED);
    }
    SUBCASE("THP never with a large anonymous set -> disabled, ahead of a process disable") {
        bwa_thp_usage_t a = usage(27600 * 1024LL, 0, 0, 0);
        CHECK(bwa_thp_assess(NV, 1, &a) == BWA_THP_WARN_DISABLED);
        CHECK(bwa_thp_assess(NV, 0, &a) == BWA_THP_WARN_DISABLED);
    }
    SUBCASE("THP_enabled 0 with a large anonymous set -> process disabled, not fallback") {
        bwa_thp_usage_t a = usage(27600 * 1024LL, 0, 0, 0);
        CHECK(bwa_thp_assess(MA, 0, &a) == BWA_THP_WARN_PROCESS_DISABLED);
        CHECK(bwa_thp_assess(AL, 0, &a) == BWA_THP_WARN_PROCESS_DISABLED);
        CHECK(bwa_thp_assess(MA, -1, &a) == BWA_THP_OK);   // unknown: eligibility decides (0 eligible)
        bwa_thp_usage_t small = usage(BWA_THP_MIN_ANON_KB - 1, 0, 0, 0);
        CHECK(bwa_thp_assess(MA, 0, &small) == BWA_THP_OK);
    }
    SUBCASE("anonymous floor: exactly BWA_THP_MIN_ANON_KB warns, 1 kB under does not") {
        bwa_thp_usage_t at = usage(BWA_THP_MIN_ANON_KB, 0);
        bwa_thp_usage_t under = usage(BWA_THP_MIN_ANON_KB - 1, 0);
        CHECK(bwa_thp_assess(NV, 1, &at) == BWA_THP_WARN_DISABLED);
        CHECK(bwa_thp_assess(MA, 1, &at) == BWA_THP_WARN_FALLBACK);
        CHECK(bwa_thp_assess(MA, 0, &at) == BWA_THP_WARN_PROCESS_DISABLED);
        CHECK(bwa_thp_assess(NV, 1, &under) == BWA_THP_OK);
        CHECK(bwa_thp_assess(MA, 1, &under) == BWA_THP_OK);
        CHECK(bwa_thp_assess(MA, 0, &under) == BWA_THP_OK);
    }
    SUBCASE("absolute floor: exactly 1 GiB short warns, 1 kB under does not") {
        bwa_thp_usage_t a = usage(5 * GiB_kB, 5 * GiB_kB - (BWA_THP_MIN_SHORT_KB - 1));
        bwa_thp_usage_t b = usage(5 * GiB_kB, 5 * GiB_kB - BWA_THP_MIN_SHORT_KB);
        CHECK(bwa_thp_assess(MA, 1, &a) == BWA_THP_OK);
        CHECK(bwa_thp_assess(MA, 1, &b) == BWA_THP_WARN_FALLBACK);
    }
    SUBCASE("relative floor: exactly 4% short warns, 1 kB under does not") {
        bwa_thp_usage_t a = usage(100 * GiB_kB, 96 * GiB_kB + 1);
        bwa_thp_usage_t b = usage(100 * GiB_kB, 96 * GiB_kB);
        CHECK(bwa_thp_assess(MA, 1, &a) == BWA_THP_OK);
        CHECK(bwa_thp_assess(MA, 1, &b) == BWA_THP_WARN_FALLBACK);
        // the relative floor is over the eligible set, not all anonymous memory
        bwa_thp_usage_t c = usage(100 * GiB_kB, 24 * GiB_kB, 26 * GiB_kB, 24 * GiB_kB);   // 2 of 26 = 7.7%
        CHECK(bwa_thp_assess(MA, 1, &c) == BWA_THP_WARN_FALLBACK);
    }
    SUBCASE("the size gate is total Anonymous, not eligible memory") {
        // 28 GiB anonymous, only 3 GiB of it eligible, 1.5 GiB of that not huge.
        bwa_thp_usage_t a = usage(28 * GiB_kB, 26 * GiB_kB + GiB_kB / 2, 3 * GiB_kB, GiB_kB + GiB_kB / 2);
        REQUIRE(a.elig_kb < BWA_THP_MIN_ANON_KB);
        CHECK(bwa_thp_assess(MA, 1, &a) == BWA_THP_WARN_FALLBACK);
        // the same eligible numbers in a process under the floor: no verdict.
        bwa_thp_usage_t b = usage(BWA_THP_MIN_ANON_KB - 1, 2 * GiB_kB, 3 * GiB_kB, GiB_kB + GiB_kB / 2);
        CHECK(bwa_thp_assess(MA, 1, &b) == BWA_THP_OK);
    }
    SUBCASE("unknown mode or NULL never warns") {
        bwa_thp_usage_t a = usage(27600 * 1024LL, 0);
        CHECK(bwa_thp_assess(BWA_THP_MODE_UNKNOWN, 1, &a) == BWA_THP_OK);
        CHECK(bwa_thp_assess(BWA_THP_MODE_INHERIT, 1, &a) == BWA_THP_OK);
        CHECK(bwa_thp_assess(MA, 1, NULL) == BWA_THP_OK);
    }
    SUBCASE("the shortfall is clamped at 0") {
        CHECK(bwa_thp_short_kb(27 * GiB_kB, 28 * GiB_kB) == 0);
        CHECK(bwa_thp_short_kb(28 * GiB_kB, 27 * GiB_kB) == GiB_kB);
    }
}

TEST_CASE("thp: sampling schedule is powers of two" * doctest::test_suite("unit/bwa_thp")) {
    CHECK_FALSE(bwa_thp_should_sample(0));
    CHECK_FALSE(bwa_thp_should_sample(-4));
    CHECK(bwa_thp_should_sample(1));
    CHECK(bwa_thp_should_sample(2));
    CHECK_FALSE(bwa_thp_should_sample(3));
    CHECK(bwa_thp_should_sample(4));
    CHECK_FALSE(bwa_thp_should_sample(6));
    CHECK(bwa_thp_should_sample(1024));
    int n = 0;
    for (long b = 1; b <= 100000; ++b) n += bwa_thp_should_sample(b);
    CHECK(n == 17);   // 2^0 .. 2^16
}

// A scripted /proc and /sys for bwa_thp_step().
struct FakeIo {
    long page = 4096;
    std::map<std::string, std::string> files;   // path -> text; missing = absent
    std::vector<bwa_thp_usage_t> walks;          // consumed in order; the last repeats
    std::vector<std::string> statuses;           // consumed with the walks; the last repeats
    size_t walk_i = 0, status_i = 0;
    int walk_fails = 0;
    int n_walks = 0, n_status = 0, n_knob_reads = 0;
    std::string buf;

    static long page_size(void *c) { return static_cast<FakeIo *>(c)->page; }
    static const char *read_text(void *c, const char *path) {
        FakeIo *f = static_cast<FakeIo *>(c);
        if (std::strcmp(path, "/proc/self/status") == 0) {
            f->n_status++;
            if (f->statuses.empty()) return NULL;
            f->buf = f->statuses[std::min(f->status_i++, f->statuses.size() - 1)];
            return f->buf.c_str();
        }
        f->n_knob_reads++;
        auto it = f->files.find(path);
        if (it == f->files.end()) return NULL;
        f->buf = it->second;
        return f->buf.c_str();
    }
    static int walk(void *c, bwa_thp_usage_t *out) {
        FakeIo *f = static_cast<FakeIo *>(c);
        f->n_walks++;
        if (f->walk_fails) return -1;
        *out = f->walks[std::min(f->walk_i++, f->walks.size() - 1)];
        return 0;
    }
    bwa_thp_io_t io() { return { this, page_size, read_text, walk }; }
};

static const char *GLOBAL = BWA_THP_KNOB_GLOBAL;
static const char *PMD = BWA_THP_KNOB_PMD;

TEST_CASE("thp: per-batch step" * doctest::test_suite("unit/bwa_thp")) {
    bwa_thp_monitor_t m;
    std::memset(&m, 0, sizeof m);
    FakeIo f;
    f.files[GLOBAL] = "always [madvise] never\n";
    const bwa_thp_usage_t healthy = usage(28 * GiB_kB, 28 * GiB_kB);
    const bwa_thp_usage_t frag = usage(28 * GiB_kB, 25 * GiB_kB);
    const std::string st_ok = status(28 * GiB_kB, 1);
    bwa_thp_io_t io = f.io();
    bwa_thp_report_t r;

    SUBCASE("samples batches 1, 2, 4, 8 and warns at most once") {
        f.walks = { healthy, healthy, frag };
        f.statuses = { st_ok };
        std::vector<int> sampled;
        int warnings = 0;
        for (int b = 1; b <= 40; ++b) {
            const bwa_thp_verdict_t v = bwa_thp_step(&m, 3, &io, &r);
            if (r.sampled) sampled.push_back(b);
            if (v != BWA_THP_OK) {
                ++warnings;
                CHECK(v == BWA_THP_WARN_FALLBACK);
                CHECK(b == 4);
            }
        }
        CHECK(warnings == 1);
        CHECK(sampled == std::vector<int>({ 1, 2, 4 }));
        CHECK(m.done == 1);
        CHECK(f.n_walks == 3);
        CHECK(m.cfg.pmd_mode == BWA_THP_MODE_MADVISE);
        CHECK(f.n_knob_reads == 9);   // config read once: global, 2048kB, 7 other sizes
    }
    SUBCASE("verbose 2 (warnings) still checks and warns") {
        f.walks = { frag };
        f.statuses = { st_ok };
        CHECK(bwa_thp_step(&m, 2, &io, &r) == BWA_THP_WARN_FALLBACK);
        CHECK(r.stop == BWA_THP_STOP_NONE);
    }
    SUBCASE("verbose < 2 stops without reading anything") {
        f.walks = { frag };
        f.statuses = { st_ok };
        CHECK(bwa_thp_step(&m, 1, &io, &r) == BWA_THP_OK);
        CHECK(m.done == 1);
        CHECK(r.sampled == 0);
        for (int b = 0; b < 8; ++b) bwa_thp_step(&m, 3, &io, &r);
        CHECK(f.n_walks == 0);
        CHECK(f.n_knob_reads == 0);
    }
    SUBCASE("a non-4 KB base page or an unknown mode stops before walking") {
        f.walks = { frag };
        f.statuses = { st_ok };
        bwa_thp_stop_t want = BWA_THP_STOP_CONFIG;
        SUBCASE("64 KB pages") { f.page = 65536; want = BWA_THP_STOP_PAGE_SIZE; }
        SUBCASE("unreadable enabled") { f.files.clear(); }
        SUBCASE("unparseable enabled") { f.files[GLOBAL] = "always madvise never\n"; }
        SUBCASE("a smaller mTHP size enabled") {
            f.files["/sys/kernel/mm/transparent_hugepage/hugepages-64kB/enabled"] =
                "always inherit [madvise] never\n";
        }
        SUBCASE("an unparseable smaller mTHP size") {
            f.files["/sys/kernel/mm/transparent_hugepage/hugepages-64kB/enabled"] = "garbage\n";
        }
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_OK);
        CHECK(r.sampled == 1);
        CHECK(r.stop == want);
        CHECK(m.done == 1);
        CHECK(f.n_walks == 0);
    }
    SUBCASE("unreadable smaps stops") {
        f.walk_fails = 1;
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_OK);
        CHECK(r.stop == BWA_THP_STOP_SMAPS);
        CHECK(m.done == 1);
        CHECK(f.n_walks == 1);
    }
    SUBCASE("stop reasons have text; NONE has none") {
        CHECK(std::string(bwa_thp_stop_reason(BWA_THP_STOP_NONE)).empty());
        CHECK(std::string(bwa_thp_stop_reason(BWA_THP_STOP_PAGE_SIZE)) == "base page is not 4 KB");
        CHECK(std::string(bwa_thp_stop_reason(BWA_THP_STOP_CONFIG)).find("THP mode") == 0);
        CHECK(std::string(bwa_thp_stop_reason(BWA_THP_STOP_SMAPS)) == "/proc/self/smaps unreadable");
    }
    SUBCASE("an incomplete walk is re-walked, up to 3 times per sample") {
        const bwa_thp_usage_t partial = usage(10 * GiB_kB, 10 * GiB_kB);
        f.walks = { partial, frag };
        f.statuses = { st_ok };
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_WARN_FALLBACK);
        CHECK(r.walks == 2);
        CHECK(r.complete == 1);
        CHECK(r.stop == BWA_THP_STOP_NONE);
    }
    SUBCASE("an incomplete sample is retried on the next batch once, then skipped") {
        const bwa_thp_usage_t partial = usage(10 * GiB_kB, 5 * GiB_kB);
        f.walks = { partial };
        f.statuses = { status(24 * GiB_kB, 1) };
        // batch 1 (scheduled): 3 incomplete walks; batch 2 is scheduled anyway.
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_OK);
        CHECK((r.sampled == 1 && r.walks == 3 && r.complete == 0));
        CHECK(m.retry_pending == 1);
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_OK);   // batch 2: scheduled
        CHECK((r.sampled == 1 && r.retry == 0));
        CHECK(m.retry_pending == 1);
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_OK);   // batch 3: the deferred retry
        CHECK((r.sampled == 1 && r.retry == 1 && r.complete == 0));
        CHECK(m.retry_pending == 0);                          // only one retry
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_OK);   // batch 4: scheduled
        CHECK(r.sampled == 1);
        CHECK(m.retry_pending == 1);
        // batch 5: the retry now sees a complete walk and warns.
        f.walks = { frag };
        f.walk_i = 0;
        f.statuses = { st_ok };
        f.status_i = 0;
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_WARN_FALLBACK);
        CHECK(r.retry == 1);
        CHECK(m.incomplete == 4);
        CHECK(m.done == 1);
        // batches 6, 7 would not have been sampled without the retry.
        bwa_thp_monitor_t m2;
        std::memset(&m2, 0, sizeof m2);
        f.walks = { partial };
        f.walk_i = 0;
        f.statuses = { status(24 * GiB_kB, 1) };
        f.status_i = 0;
        std::vector<int> sampled;
        for (int b = 1; b <= 9; ++b) {
            bwa_thp_step(&m2, 3, &io, &r);
            if (r.sampled) sampled.push_back(b);
        }
        CHECK(sampled == std::vector<int>({ 1, 2, 3, 4, 5, 8, 9 }));
    }
    SUBCASE("THP_enabled: 0 gives the process-disabled verdict, not fallback") {
        f.walks = { usage(28 * GiB_kB, 0, 0, 0) };
        f.statuses = { status(28 * GiB_kB, 0) };
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_WARN_PROCESS_DISABLED);
        CHECK(m.done == 1);
    }
    SUBCASE("unadvised memory in madvise mode does not warn") {
        f.walks = { parse_smaps(vma("7f0000000000-7f0400000000", 10 * GiB_kB, 10 * GiB_kB, 1)
                                + vma("7f0400000000-7f0800000000", 16 * GiB_kB, 0, 0)) };
        f.statuses = { status(26 * GiB_kB, 1) };
        for (int b = 1; b <= 16; ++b) CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_OK);
        CHECK(m.done == 0);
        CHECK(f.n_walks == 5);
    }
    SUBCASE("global never with no per-size knob -> disabled naming the global knob") {
        f.files[GLOBAL] = "always madvise [never]\n";
        f.walks = { usage(28 * GiB_kB, 0, 0, 0) };
        f.statuses = { st_ok };
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_WARN_DISABLED);
        CHECK(std::string(m.cfg.knob) == GLOBAL);
    }
    SUBCASE("2 MB never via hugepages-2048kB -> disabled naming that knob") {
        f.files[PMD] = "always inherit madvise [never]\n";
        f.walks = { usage(28 * GiB_kB, 0, 0, 0) };
        f.statuses = { st_ok };
        CHECK(bwa_thp_step(&m, 3, &io, &r) == BWA_THP_WARN_DISABLED);
        CHECK(std::string(m.cfg.knob) == PMD);
    }
    SUBCASE("NULL monitor or io is a no-op") {
        CHECK(bwa_thp_step(NULL, 3, &io, &r) == BWA_THP_OK);
        CHECK(bwa_thp_step(&m, 3, NULL, &r) == BWA_THP_OK);
        CHECK(r.sampled == 0);
    }
}

static int count_lines(const std::string &s)
{
    int n = 0;
    for (char c : s) n += c == '\n';
    return n;
}

TEST_CASE("thp: warning text" * doctest::test_suite("unit/bwa_thp")) {
    char buf[1024];
    const bwa_thp_verdict_t verdicts[3] = { BWA_THP_WARN_FALLBACK, BWA_THP_WARN_DISABLED,
                                            BWA_THP_WARN_PROCESS_DISABLED };
    SUBCASE("every verdict: two '[W::f] ' lines, return value == strlen") {
        bwa_thp_usage_t u = usage(28 * GiB_kB, 25 * GiB_kB);
        for (bwa_thp_verdict_t v : verdicts) {
            const int r = bwa_thp_format_warning(v, &u, NULL, "f", buf, sizeof buf);
            REQUIRE(r > 0);
            REQUIRE(r < (int)sizeof buf);
            CHECK(r == (int)std::strlen(buf));
            const std::string s(buf);
            CHECK(count_lines(s) == 2);
            CHECK(s.back() == '\n');
            CHECK(s.rfind("[W::f] ", 0) == 0);
            CHECK(s.substr(s.find('\n') + 1).rfind("[W::f] ", 0) == 0);
            CHECK(s.find("on human WGS this costs on the order of 1% wall time per GiB")
                  != std::string::npos);
            CHECK(s.find("see docs/src/best-practices/optimization-checklist.md, "
                         "\"Keep transparent huge pages working\".\n")
                  != std::string::npos);
        }
    }
    SUBCASE("fallback reports the eligible shortfall and a typical cause") {
        bwa_thp_usage_t u = usage(30 * GiB_kB, 25 * GiB_kB, 28 * GiB_kB, 25 * GiB_kB);
        bwa_thp_format_warning(BWA_THP_WARN_FALLBACK, &u, NULL, "f", buf, sizeof buf);
        const std::string s(buf);
        CHECK(s.rfind("[W::f] 3.0 GiB of 28.0 GiB of THP-eligible anonymous memory is on 4 KB pages", 0) == 0);
        CHECK(s.find("typically because host memory was full of page cache or fragmented")
              != std::string::npos);
        CHECK(s.find("\n[W::f] start with free memory") != std::string::npos);
    }
    SUBCASE("disabled reports the clamped shortfall, not total anonymous, and names the knob") {
        bwa_thp_usage_t u = usage(28 * GiB_kB, 10 * GiB_kB);
        bwa_thp_format_warning(BWA_THP_WARN_DISABLED, &u, NULL, "f", buf, sizeof buf);
        std::string s(buf);
        CHECK(s.find("(/sys/kernel/mm/transparent_hugepage/enabled is [never]): 18.0 GiB of 28.0 GiB")
              != std::string::npos);
        CHECK(s.find("(about 28% on 5M-pair human WGS)") != std::string::npos);
        CHECK(s.find("\n[W::f] set it to madvise or always") != std::string::npos);
        bwa_thp_format_warning(BWA_THP_WARN_DISABLED, &u, PMD, "f", buf, sizeof buf);
        s = buf;
        CHECK(s.find(std::string("(") + PMD + " is [never])") != std::string::npos);
    }
    SUBCASE("process-disabled names the prctl cause and the remedy") {
        bwa_thp_usage_t u = usage(27 * GiB_kB, 0);
        bwa_thp_format_warning(BWA_THP_WARN_PROCESS_DISABLED, &u, NULL, "f", buf, sizeof buf);
        const std::string s(buf);
        CHECK(s.rfind("[W::f] transparent huge pages are disabled for this process", 0) == 0);
        CHECK(s.find("THP_enabled: 0") != std::string::npos);
        CHECK(s.find("MIMALLOC_ALLOW_THP=0") != std::string::npos);
        CHECK(s.find("27.0 GiB of 27.0 GiB") != std::string::npos);
        CHECK(s.find("THP-eligible") == std::string::npos);   // not the FALLBACK text
    }
    SUBCASE("AnonHugePages > Anonymous prints 0.0 GiB, not a negative") {
        bwa_thp_usage_t u = usage(27 * GiB_kB, 28 * GiB_kB);
        bwa_thp_format_warning(BWA_THP_WARN_FALLBACK, &u, NULL, "f", buf, sizeof buf);
        CHECK(std::string(buf).find("0.0 GiB of 27.0 GiB") != std::string::npos);
        bwa_thp_format_warning(BWA_THP_WARN_DISABLED, &u, NULL, "f", buf, sizeof buf);
        CHECK(std::string(buf).find("]): 0.0 GiB of 27.0 GiB") != std::string::npos);
    }
    SUBCASE("NULL func defaults to mem") {
        bwa_thp_usage_t u = usage(28 * GiB_kB, 25 * GiB_kB);
        bwa_thp_format_warning(BWA_THP_WARN_FALLBACK, &u, NULL, NULL, buf, sizeof buf);
        const std::string s(buf);
        CHECK(s.rfind("[W::mem] ", 0) == 0);
        CHECK(s.find("\n[W::mem] ") != std::string::npos);
    }
    SUBCASE("a small buffer truncates and reports the full length") {
        bwa_thp_usage_t u = usage(28 * GiB_kB, 25 * GiB_kB);
        const int full = bwa_thp_format_warning(BWA_THP_WARN_FALLBACK, &u, NULL, "f", buf, sizeof buf);
        char small[16];
        const int r = bwa_thp_format_warning(BWA_THP_WARN_FALLBACK, &u, NULL, "f", small, sizeof small);
        CHECK(r == full);
        CHECK(std::strlen(small) == sizeof small - 1);
        CHECK(std::string(small) == std::string(buf).substr(0, sizeof small - 1));
        CHECK(bwa_thp_format_warning(BWA_THP_WARN_FALLBACK, &u, NULL, "f", small, 0) == 0);
        CHECK(bwa_thp_format_warning(BWA_THP_WARN_FALLBACK, &u, NULL, "f", NULL, 8) == 0);
    }
    SUBCASE("NULL usage or OK writes nothing") {
        bwa_thp_usage_t u = usage(1, 1);
        buf[0] = 'x';
        CHECK(bwa_thp_format_warning(BWA_THP_WARN_FALLBACK, NULL, NULL, "f", buf, sizeof buf) == 0);
        CHECK(buf[0] == '\0');
        buf[0] = 'x';
        CHECK(bwa_thp_format_warning(BWA_THP_OK, &u, NULL, "f", buf, sizeof buf) == 0);
        CHECK(buf[0] == '\0');
    }
}
