// Unit tests for bns_meth_seed_mismatch: a --meth seed index (<prefix>.meth) must be the converted
// copy `index --meth` writes for the original reference (<prefix>), contigs r<name>, f<name> per
// original contig, in order, each the original's length. A plain `index` on the same prefix
// rewrites only the original set; mem --meth must refuse the mix rather than map seeds onto
// whatever original contig sits at rid / 2.

#include "doctest/doctest.h"

#include "bntseq.h"

#include <string>
#include <vector>

namespace {

struct Bns {
    std::vector<std::string> names;
    std::vector<bntann1_t> anns;
    bntseq_t b{};
    Bns(const std::vector<std::pair<std::string, int32_t>> &contigs)
    {
        names.reserve(contigs.size());
        int64_t off = 0;
        for (const auto &c : contigs) {
            names.push_back(c.first);
            bntann1_t a{};
            a.offset = off;
            a.len = c.second;
            off += c.second;
            anns.push_back(a);
        }
        for (size_t i = 0; i < anns.size(); ++i) anns[i].name = &names[i][0];
        b.l_pac = off;
        b.n_seqs = (int32_t)anns.size();
        b.anns = anns.data();
    }
    Bns(const Bns &) = delete;
    Bns &operator=(const Bns &) = delete;
};

std::vector<std::pair<std::string, int32_t>> seed_of(const std::vector<std::pair<std::string, int32_t>> &o)
{
    std::vector<std::pair<std::string, int32_t>> s;
    for (const auto &c : o) {
        s.push_back({"r" + c.first, c.second});
        s.push_back({"f" + c.first, c.second});
    }
    return s;
}

} // namespace

TEST_SUITE_BEGIN("unit/bns_meth_seed_mismatch");

TEST_CASE("bns_meth_seed_mismatch: the seed index index --meth writes matches its original")
{
    const std::vector<std::pair<std::string, int32_t>> orig = {{"chr1", 5000}, {"chr2", 300}, {"chrM", 16}};
    Bns o(orig), s(seed_of(orig));
    CHECK(bns_meth_seed_mismatch(&s.b, &o.b) == nullptr);
}

TEST_CASE("bns_meth_seed_mismatch: a seed index from another build is refused")
{
    const std::vector<std::pair<std::string, int32_t>> orig = {{"chr1", 5000}, {"chr2", 300}};
    Bns o(orig);
    {   // one contig fewer
        Bns s(seed_of({{"chr1", 5000}}));
        CHECK(std::string(bns_meth_seed_mismatch(&s.b, &o.b)) == "contig count");
    }
    {   // same count and total, one contig resized
        auto x = seed_of(orig);
        x[0].second -= 1;
        x[2].second += 1;
        Bns s(x);
        CHECK(std::string(bns_meth_seed_mismatch(&s.b, &o.b)) == "contig length");
    }
    {   // a longer reference
        auto x = seed_of(orig);
        x[3].second += 7;
        Bns s(x);
        CHECK(std::string(bns_meth_seed_mismatch(&s.b, &o.b)) == "reference length");
    }
    {   // same lengths, renamed contig
        Bns s(seed_of({{"chr1", 5000}, {"chrX", 300}}));
        CHECK(std::string(bns_meth_seed_mismatch(&s.b, &o.b)) == "contig name");
    }
    {   // the two converted copies swapped (f before r)
        auto x = seed_of(orig);
        std::swap(x[0], x[1]);
        Bns s(x);
        CHECK(std::string(bns_meth_seed_mismatch(&s.b, &o.b)) == "contig name");
    }
}

TEST_CASE("meth metadata checks the forward copy as well as the reverse")
{
    const std::vector<std::pair<std::string, int32_t>> orig = {{"chr1", 5000}, {"chr2", 300}};
    Bns o(orig);
    auto names = seed_of(orig);
    names[1].first = "fwrong";
    Bns bad_name(names);
    CHECK(std::string(bns_meth_seed_mismatch(&bad_name.b, &o.b)) == "contig name");
    auto lengths = seed_of(orig);
    lengths[1].second--;
    lengths[3].second++;
    Bns bad_length(lengths);
    CHECK(std::string(bns_meth_seed_mismatch(&bad_length.b, &o.b)) == "contig length");
    o.b.l_pac = INT64_MAX;
    CHECK(std::string(bns_meth_seed_mismatch(&bad_length.b, &o.b)) == "reference length");
}

TEST_CASE("meth metadata accepts upstream empty and truncated name round trips")
{
    Bns empty_original({{"(null)", 5000}}), empty_seed({{"r", 5000}, {"f", 5000}});
    CHECK(bns_meth_seed_mismatch(&empty_seed.b, &empty_original.b) == nullptr);
    const std::string long_name(8192, 'x');
    Bns long_original({{long_name, 5000}});
    Bns long_seed({{"r" + long_name.substr(0, 8191), 5000},
                   {"f" + long_name.substr(0, 8191), 5000}});
    CHECK(bns_meth_seed_mismatch(&long_seed.b, &long_original.b) == nullptr);
}

TEST_SUITE_END();
