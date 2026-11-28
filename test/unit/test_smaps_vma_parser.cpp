// Unit tests for the streaming one-VMA smaps parser (bwa_smaps_vma_*, src/bwa_thp.cpp) that the shm
// huge-page check reads /proc/self/smaps through, in chunks of any size. Parser inputs are canned
// smaps text in the shape of a 6.1 kernel's. The verdict must not depend on where the chunk
// boundaries fall (lines span chunks), must read only the VMA holding the address, and must
// agree with the THP parser on what a VMA header is.

#include "doctest/doctest.h"

#include "bwa_thp.h"
#include "temp_dir.h"

#include <algorithm>
#include <cstdlib>
#include <string>

namespace {

// One VMA record. shmem_pmd_kb / file_pmd_kb < 0 leave that field out (an old kernel prints
// neither). A hugetlb mapping has kernel_page_kb 2048 or more and its pages in Shared_Hugetlb,
// with Rss 0.
std::string vma(const char *range, long long rss_kb, long long shmem_pmd_kb, long long file_pmd_kb,
                long long kernel_page_kb = 4, long long shared_hugetlb_kb = 0)
{
    std::string s = std::string(range) + " r--s 00000000 00:19 1234     /dev/shm/bwa_idx_a_very_long_name\n";
    s += "Size:           10485760 kB\n"
         "KernelPageSize: " + std::to_string(kernel_page_kb) + " kB\n"
         "MMUPageSize:    " + std::to_string(kernel_page_kb) + " kB\n"
         "Rss:            " + std::to_string(rss_kb) + " kB\n"
         "Shared_Clean:          0 kB\n"
         "AnonHugePages:         0 kB\n";
    if (shmem_pmd_kb >= 0) s += "ShmemPmdMapped: " + std::to_string(shmem_pmd_kb) + " kB\n";
    if (file_pmd_kb >= 0) s += "FilePmdMapped:  " + std::to_string(file_pmd_kb) + " kB\n";
    s += "Shared_Hugetlb: " + std::to_string(shared_hugetlb_kb) + " kB\n"
         "Private_Hugetlb:       0 kB\n"
         "THPeligible:    0\n"
         "VmFlags: rd sh mr mw me ms hg \n";
    return s;
}

struct Verdict { int v; int done; };

Verdict chunked(const std::string &text, unsigned long long addr, size_t chunk)
{
    bwa_smaps_vma_parser_t p;
    bwa_smaps_vma_init(&p, addr);
    int done = 0;
    for (size_t i = 0; i < text.size() && !done; i += chunk) {
        bwa_smaps_vma_feed(&p, text.data() + i, std::min(chunk, text.size() - i));
        done = bwa_smaps_vma_done(&p);
    }
    return Verdict{ bwa_smaps_vma_finish(&p), done };
}

const unsigned long long SEG = 0x7f0000000000ULL;
const unsigned long long NEXT = 0x7f0300000000ULL;   // the VMA after SEG, on 2 MB pages

std::string around(const std::string &target)
{
    return vma("55d000000000-55d000021000", 132, 0, 0) + target
         + vma("7f0300000000-7f0300200000", 2048, 2048, 0);
}

#if defined(__linux__)
const bool proc_reader_unavailable = false;
#else
const bool proc_reader_unavailable = true;
#endif

// Keep the allocation owned even if a REQUIRE exits the test early.
struct TextBuffer {
    char *data = NULL;
    size_t cap = 0;
    ~TextBuffer() { std::free(data); }
};

} // namespace

TEST_CASE("smaps vma parser: the verdict does not depend on chunking"
          * doctest::test_suite("unit/smaps_vma_parser"))
{
    struct Case { long long rss, shmem, file, kpage, hugetlb; int want; } cases[] = {
        {2048, 2048, 0, 4, 0, 1},       // ShmemPmdMapped
        {2048, 0, 2048, 4, 0, 1},       // FilePmdMapped
        {4096, 0, 0, 4, 0, 0},          // resident, 4 KB pages
        {0, 0, 0, 4, 0, -1},            // nothing resident: cannot tell
        {4096, -1, -1, 4, 0, -1},       // a kernel that prints neither PMD field
        {0, 0, 0, 1048576, 10485760, 1},   // hugetlb, 1 GB pages: Rss 0, pages in Shared_Hugetlb
        {0, 0, 0, 4, 10485760, 1},      // hugetlb pages, whatever KernelPageSize says
        {0, 0, 0, 2048, 0, 1},          // a hugetlb mapping not yet faulted in
    };
    for (const Case &c : cases) {
        CAPTURE(c.want);
        CAPTURE(c.kpage);
        const std::string text = around(vma("7f0000000000-7f0280000000", c.rss, c.shmem, c.file,
                                            c.kpage, c.hugetlb));
        for (size_t chunk : {(size_t)1, (size_t)2, (size_t)7, (size_t)63, (size_t)64, (size_t)65, text.size()}) {
            CAPTURE(chunk);
            CHECK(chunked(text, SEG, chunk).v == c.want);
            CHECK(chunked(text, SEG + 0x1000, chunk).v == c.want);   // any address inside it
        }
        // An address no VMA holds (the end is exclusive), and one in the next VMA (which is huge).
        CHECK(chunked(text, 0x10ULL, 5).v == -1);
        CHECK(chunked(text, 0x7f0280000000ULL, 5).v == -1);
        CHECK(chunked(text, NEXT, 5).v == 1);
        // A last line without its newline is still read.
        CHECK(chunked(text.substr(0, text.size() - 1), NEXT, 3).v == 1);
    }
    CHECK(chunked("", SEG, 1).v == -1);
}

TEST_CASE("smaps vma parser: the neighbours' fields are not read, and the read can stop after the VMA"
          * doctest::test_suite("unit/smaps_vma_parser"))
{
    // Ours is on 4 KB pages; the VMA after it is PMD-mapped.
    const std::string text = around(vma("7f0000000000-7f0280000000", 4096, 0, 0));
    const Verdict ours = chunked(text, SEG, 16);
    CHECK(ours.v == 0);
    CHECK(ours.done == 1);   // stopped at the next VMA's header
    const Verdict last = chunked(text, NEXT, 16);
    CHECK(last.v == 1);
    CHECK(last.done == 0);   // the last VMA: read to the end
}

TEST_CASE("smaps vma parser: only '<hex>-<hex>' opens a VMA"
          * doctest::test_suite("unit/smaps_vma_parser"))
{
    // Lines inside the VMA that start like a header but are not one: a hex word with no '-', a
    // bare number, "0x...". The old rule (any line starting with a digit or a-f) took each of
    // them for a header, which would close the VMA before its FilePmdMapped line.
    const std::string text = "7f0000000000-7f0280000000 r--s 00000000 00:19 1234 /dev/shm/x\n"
                             "Rss:            4096 kB\n"
                             "deadbeef 0 kB\n"
                             "7f00 r--s\n"
                             "0 kB\n"
                             "0x7f00-0x7f01\n"
                             "FilePmdMapped:  2048 kB\n";
    const Verdict v = chunked(text, SEG, 9);
    CHECK(v.v == 1);
    CHECK(v.done == 0);
    bwa_thp_smaps_parser_t t;
    bwa_thp_smaps_init(&t);
    bwa_thp_smaps_feed(&t, text.data(), text.size());
    bwa_thp_usage_t u;
    bwa_thp_smaps_finish(&t, &u);
    CHECK(u.vmas == 1);
}

// These exercise the real file reader on Linux; other hosts report explicit skips.
TEST_CASE("proc text reader: a smaller limit applies to retained capacity"
          * doctest::test_suite("unit/smaps_vma_parser")
          * doctest::skip(proc_reader_unavailable))
{
    const bwa_tests::TempDir dir("proc_text");
    const std::string path = dir.file("status");
    TextBuffer b;
    dir.write("status", std::string(4095, 'a'));
    REQUIRE(bwa_proc_read_text(path.c_str(), &b.data, &b.cap, 4096) == 4095);
    char *const retained = b.data;
    const size_t capacity = b.cap;
    REQUIRE(capacity >= 4096);

    const std::string text(100, 'b');
    dir.write("status", text);
    CHECK(bwa_proc_read_text(path.c_str(), &b.data, &b.cap, 64) == -1);
    CHECK(b.data == retained);
    CHECK(b.cap == capacity);
    REQUIRE(bwa_proc_read_text(path.c_str(), &b.data, &b.cap, 101) == 100);
    CHECK(b.data == retained);
    CHECK(b.cap == capacity);
    CHECK(std::string(b.data) == text);
    CHECK(b.data[100] == '\0');
}

TEST_CASE("proc text reader: the bound includes the terminating NUL"
          * doctest::test_suite("unit/smaps_vma_parser")
          * doctest::skip(proc_reader_unavailable))
{
    const bwa_tests::TempDir dir("proc_text");
    const std::string path = dir.file("status");
    TextBuffer b;
    dir.write("status", std::string(4095, 'a'));
    REQUIRE(bwa_proc_read_text(path.c_str(), &b.data, &b.cap, 4096) == 4095);
    dir.write("status", std::string(63, 'c'));
    REQUIRE(bwa_proc_read_text(path.c_str(), &b.data, &b.cap, 64) == 63);
    CHECK(std::string(b.data) == std::string(63, 'c'));
    CHECK(b.data[63] == '\0');
    dir.write("status", std::string(64, 'c'));
    CHECK(bwa_proc_read_text(path.c_str(), &b.data, &b.cap, 64) == -1);
}

TEST_CASE("proc text reader: an empty file still needs one byte"
          * doctest::test_suite("unit/smaps_vma_parser")
          * doctest::skip(proc_reader_unavailable))
{
    const bwa_tests::TempDir dir("proc_text");
    const std::string path = dir.file("status");
    TextBuffer b;
    dir.write("status", "");
    CHECK(bwa_proc_read_text(path.c_str(), &b.data, &b.cap, 0) == -1);
    REQUIRE(bwa_proc_read_text(path.c_str(), &b.data, &b.cap, 1) == 0);
    CHECK(b.data[0] == '\0');
    CHECK(bwa_proc_read_text(path.c_str(), &b.data, &b.cap, 0) == -1);
}
