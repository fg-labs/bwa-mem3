// Keep htslib's kstring type out of translation units that use the aligner's type.
#include "htslib/sam.h"
#include "doctest/doctest.h"
#include "aln_bam_decode.h"

#include <cstdlib>
#include <memory>

namespace bwa_tests {
std::string decoded_aln_records(const std::string &path)
{
    std::unique_ptr<htsFile, decltype(&hts_close)> input(hts_open(path.c_str(), "r"), hts_close);
    REQUIRE(input != nullptr);
    std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> header(sam_hdr_read(input.get()), sam_hdr_destroy);
    REQUIRE(header != nullptr);
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> record(bam_init1(), bam_destroy1);
    REQUIRE(record != nullptr);
    std::string records;
    int status;
    while ((status = sam_read1(input.get(), header.get(), record.get())) >= 0) {
        kstring_t text = {0, 0, nullptr};
        const int formatted = sam_format1(header.get(), record.get(), &text);
        std::unique_ptr<char, decltype(&free)> buffer(text.s, free);
        REQUIRE(formatted >= 0);
        records.append(text.s, text.l);
        records += "\n";
    }
    CHECK(status == -1); // Distinguish EOF from a truncated or malformed BAM.
    CHECK(hts_close(input.release()) == 0);
    return records;
}
}
