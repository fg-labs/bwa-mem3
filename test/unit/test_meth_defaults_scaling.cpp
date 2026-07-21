// Unit tests for mem_opt_apply_meth_defaults — the --meth scoring defaults.
//
// bwa scales score-derived options by the match score: `update_a()` multiplies
// T, pen_clip5/3, pen_unpaired, b, gap penalties and zdrop by opt->a whenever
// the user passes -A and has not set the option explicitly. Every one of those
// values is expressed in units of the match score, so a bare constant is only
// meaningful at a == 1.
//
// The --meth defaults are borrowed verbatim from bwameth
// (`bwa-mem2 mem -T 40 -B 2 -L 10 -CM`, plus -U 100 for paired), and bwameth
// runs bwa at its default a == 1. Applying those constants FLAT therefore
// silently discards -A: `--meth -A 2` produced T=40 where every other
// score-derived option had doubled, making the output threshold half as strict
// relative to the alignment scores it is compared against.
//
// Contract pinned here:
//   * a == 1  -> exactly bwameth's constants (byte-identical to prior behaviour)
//   * a  > 1  -> the same constants scaled by a (self-consistent scoring)
//   * user-supplied values (opt0 sentinel set) always win, at any a

#include "doctest/doctest.h"
#include "bwamem.h"

#include <cstdlib>

namespace {

// Fresh default opt with the match score forced to `a`, mirroring how main_mem
// has already parsed -A by the time the meth defaults are applied.
mem_opt_t *opt_with_a(int a)
{
    mem_opt_t *o = mem_opt_init();
    o->a = a;
    o->meth_mode = 1;
    return o;
}

}  // namespace

TEST_CASE("meth defaults at a=1 are bwameth's constants verbatim") {
    mem_opt_t *o = opt_with_a(1);
    mem_opt_t o0;                      // no sentinel set => nothing user-supplied
    memset(&o0, 0, sizeof(o0));
    o->meth_scoring = MEM_METH_SCORING_COLLAPSED;

    mem_opt_apply_meth_defaults(o, &o0);

    CHECK(o->T            == 40);      // bwameth -T 40
    CHECK(o->pen_clip5    == 10);      // bwameth -L 10
    CHECK(o->pen_clip3    == 10);
    CHECK(o->pen_unpaired == 100);     // bwameth -U 100
    CHECK(o->b            == 2);       // bwameth -B 2 (collapsed only)
    free(o);
}

TEST_CASE("meth defaults scale with -A so scoring stays self-consistent") {
    mem_opt_t *o = opt_with_a(2);
    mem_opt_t o0;
    memset(&o0, 0, sizeof(o0));
    o->meth_scoring = MEM_METH_SCORING_COLLAPSED;

    mem_opt_apply_meth_defaults(o, &o0);

    CHECK(o->T            == 80);      // 40 * 2
    CHECK(o->pen_clip5    == 20);      // 10 * 2
    CHECK(o->pen_clip3    == 20);
    CHECK(o->pen_unpaired == 200);     // 100 * 2
    CHECK(o->b            == 4);       // 2 * 2
    free(o);
}

TEST_CASE("explicit user values are never overwritten, at any match score") {
    mem_opt_t *o = opt_with_a(3);
    mem_opt_t o0;
    memset(&o0, 0, sizeof(o0));
    // User passed -T/-L/-U/-B explicitly: sentinels set, values already parsed.
    o0.T = 1;            o->T            = 55;
    o0.pen_clip5 = 1;    o->pen_clip5    = 7;
    o0.pen_clip3 = 1;    o->pen_clip3    = 8;
    o0.pen_unpaired = 1; o->pen_unpaired = 9;
    o0.b = 1;            o->b            = 6;
    o->meth_scoring = MEM_METH_SCORING_COLLAPSED;

    mem_opt_apply_meth_defaults(o, &o0);

    CHECK(o->T            == 55);
    CHECK(o->pen_clip5    == 7);
    CHECK(o->pen_clip3    == 8);
    CHECK(o->pen_unpaired == 9);
    CHECK(o->b            == 6);
    free(o);
}

TEST_CASE("genomic scoring keeps bwa's variant-aware mismatch, scaled") {
    // COLLAPSED adopts bwameth's lenient -B 2; GENOMIC deliberately keeps bwa's
    // default b (4 at a=1), which update_a has already scaled. So the meth
    // block must not touch b under GENOMIC.
    mem_opt_t *o = opt_with_a(2);
    mem_opt_t o0;
    memset(&o0, 0, sizeof(o0));
    o->meth_scoring = MEM_METH_SCORING_GENOMIC;
    o->b = 8;                          // as update_a would leave it (4 * 2)

    mem_opt_apply_meth_defaults(o, &o0);

    CHECK(o->b == 8);                  // untouched by the meth block
    CHECK(o->T == 80);                 // other defaults still scale
    free(o);
}

TEST_CASE("meth defaults set -M (no multi) regardless of match score") {
    mem_opt_t *o = opt_with_a(2);
    mem_opt_t o0;
    memset(&o0, 0, sizeof(o0));
    o->meth_scoring = MEM_METH_SCORING_COLLAPSED;

    mem_opt_apply_meth_defaults(o, &o0);

    CHECK((o->flag & MEM_F_NO_MULTI) != 0);
    free(o);
}
