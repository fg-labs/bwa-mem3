// Fast default, cap, gate and compatibility checks; exhaustive matrices live in integration.
#include <type_traits>

#include "bsw_compact16_test.h"

#if HAVE_BSW_VECTOR_8_16

// Case / 5 selects NEON, AVX2 or AVX-512BW; Case % 5 selects one boundary check.
TEST_CASE_TEMPLATE("lane compaction 16-bit: default and boundary checks"
                       * doctest::test_suite("unit/bandedswa"), Case,
                   std::integral_constant<int, 0>,
                   std::integral_constant<int, 1>,
                   std::integral_constant<int, 2>,
                   std::integral_constant<int, 3>,
                   std::integral_constant<int, 4>,
                   std::integral_constant<int, 5>,
                   std::integral_constant<int, 6>,
                   std::integral_constant<int, 7>,
                   std::integral_constant<int, 8>,
                   std::integral_constant<int, 9>,
                   std::integral_constant<int, 10>,
                   std::integral_constant<int, 11>,
                   std::integral_constant<int, 12>,
                   std::integral_constant<int, 13>,
                   std::integral_constant<int, 14>) {
    const int lanes[] = {8, 16, 32};
    const int selected_lanes = lanes[Case::value / 5];
    const int selected_check = Case::value % 5;
    for (const Tier &t : compacting_tiers()) {
        if (t.lanes != selected_lanes) continue;
        SUBCASE(t.name.c_str()) {
            INFO("tier " << t.name);
            std::mt19937_64 rng(0xDEFA1716ull);
            const ExtScoring sc = scoring(1, 4, 6, 1, 6, 1, 100, 5);
            const std::vector<ExtPair> batch = staggered_batch(rng, t.lanes, 4);
            std::vector<const ExtPair *> ptrs;
            for (const ExtPair &p : batch) ptrs.push_back(&p);
            uint64_t moves = 0, retired = 0;
            if (selected_check == 0) SUBCASE("the default uses the tier's configured group count") {
                const int lane_groups = BSW_COMPACT_GROUPS_MAX + 1;
                const std::vector<ExtPair> default_batch = staggered_batch(rng, t.lanes, lane_groups);
                std::vector<const ExtPair *> default_ptrs;
                for (const ExtPair &p : default_batch) default_ptrs.push_back(&p);
                std::unique_ptr<IBandedPairWiseSW> k(t.make(sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop, sc.pen_clip,
                                                            sc.mat, (int8_t)sc.a, (int8_t)sc.b, 1));
                const std::unique_ptr<IBandedPairWiseSW> plain = make(t, sc, 0, 0);
                const std::vector<ExtResult> got = bwa_tests::score_bsw_batch(*k, default_ptrs, 16, 100);
                const std::vector<ExtResult> ref = bwa_tests::score_bsw_batch(*plain, default_ptrs, 16, 100);
                CHECK(got == ref);
                const int groups = std::min(t.default_groups, BSW_COMPACT_GROUPS_MAX);
                const uint64_t expected = groups > 0 ? (uint64_t)((lane_groups + groups - 1) / groups) : 0;
                CHECK(k->lane_compaction16_counts().superblocks == expected);
                CHECK(plain->lane_compaction16_counts().superblocks == 0);
            }
            if (selected_check == 1) SUBCASE("one lane group runs the plain wrapper") {
                const std::unique_ptr<IBandedPairWiseSW> k = make(t, sc, 8, 8);
                const std::vector<const ExtPair *> one(ptrs.begin(), ptrs.begin() + t.lanes);
                bwa_tests::score_bsw_batch(*k, one, 16, 100);
                CHECK(k->lane_compaction16_counts().superblocks == 0);
            }
            if (selected_check == 2) SUBCASE("groups past the cap use the capped number of groups per superblock") {
                const std::vector<ExtPair> big = staggered_batch(rng, t.lanes, BSW_COMPACT_GROUPS_MAX + 1);
                std::vector<const ExtPair *> big_ptrs;
                for (const ExtPair &p : big) big_ptrs.push_back(&p);
                const std::unique_ptr<IBandedPairWiseSW> capped = make(t, sc, 1000, 1);
                const std::unique_ptr<IBandedPairWiseSW> plain = make(t, sc, 0, 0);
                CHECK(bwa_tests::score_bsw_batch(*capped, big_ptrs, 16, 100) ==
                      bwa_tests::score_bsw_batch(*plain, big_ptrs, 16, 100));
                CHECK(capped->lane_compaction16_counts().superblocks == 2);
            }
            if (selected_check == 3) SUBCASE("one configured group per superblock runs four superblocks without moves") {
                const std::unique_ptr<IBandedPairWiseSW> one = make(t, sc, 1, 1);
                const std::unique_ptr<IBandedPairWiseSW> plain = make(t, sc, 0, 0);
                CHECK(bwa_tests::score_bsw_batch(*one, ptrs, 16, 100) ==
                      bwa_tests::score_bsw_batch(*plain, ptrs, 16, 100));
                CHECK(one->lane_compaction16_counts().superblocks == 4);
                CHECK(one->lane_compaction16_counts().moves == 0);
            }
            if (selected_check == 4) SUBCASE("nonpositive row settings preserve result records") {
                CHECK(check_batch(t, sc, batch, 100, 1000, 0, "zero rows", moves, retired) == 0);
                CHECK(check_batch(t, sc, batch, 100, 1, -5, "negative rows", moves, retired) == 0);
            }
        }
    }
}

TEST_CASE("lane compaction 16-bit: tiers without it ignore the setting"
          * doctest::test_suite("unit/bandedswa")) {
#if defined(__x86_64__) || defined(__i386__)
    __builtin_cpu_init();
    struct Plain { const char *name; Factory make; bool ok; };
    const Plain tiers[] = {{"sse41", make_bsw_kernel_sse41, __builtin_cpu_supports("sse4.1") != 0},
                           {"sse42", make_bsw_kernel_sse42, __builtin_cpu_supports("sse4.2") != 0},
                           {"avx", make_bsw_kernel_avx, __builtin_cpu_supports("avx") != 0}};
    std::mt19937_64 rng(0x0FF16C0Dull);
    const ExtScoring sc = scoring(1, 4, 6, 1, 6, 1, 100, 5);
    for (const Plain &p : tiers) {
        if (!p.ok) continue;
        INFO("tier " << p.name);
        const Tier t{p.name, p.make, 8, 0};
        const std::vector<ExtPair> batch = staggered_batch(rng, 8, 4);
        uint64_t moves = 0, retired = 0;
        CHECK(check_batch(t, sc, batch, 100, 8, 8, "ignored", moves, retired) == 0);
        const std::unique_ptr<IBandedPairWiseSW> k = make(t, sc, 8, 8);
        std::vector<const ExtPair *> ptrs;
        for (const ExtPair &q : batch) ptrs.push_back(&q);
        (void)bwa_tests::score_bsw_batch(*k, ptrs, 16, 100);
        CHECK(k->lane_compaction16_counts().superblocks == 0);   // the compaction driver never ran
        CHECK(moves == 0);
        CHECK(retired == 0);
    }
#else
    MESSAGE("lane compaction 16-bit: every aarch64 build has the driver; nothing to check here");
#endif
}

#endif   // HAVE_BSW_VECTOR_8_16
