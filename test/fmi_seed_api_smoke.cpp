#include "fmi_seed_api.h"
#include <cstdio>
int main(int argc, char **argv) {
    CP_OCC blk; CP_OCC_PLANES pl; SMEM s; s.k = 0; s.s = 0; (void)blk; (void)pl; (void)s;
    unsigned long m = one_hot_mask_array[3];
    int pc = _mm_countbits_64(m);
    printf("CP_SHIFT=%d CP_MASK=%d popcount=%d\n", CP_SHIFT, CP_MASK, pc);

    // link-check the facade symbols exist in libbwa.a:
    FmiSeed       *(*open)(const char*)              = &fmi_seed_open;
    void           (*close)(FmiSeed*)                = &fmi_seed_close;
    const CP_OCC  *(*cp_occ)(const FmiSeed*)         = &fmi_seed_cp_occ;
    const CP_OCC_PLANES *(*cp_planes)(const FmiSeed*) = &fmi_seed_cp_occ_planes;
    const int64_t *(*cp_base)(const FmiSeed*)        = &fmi_seed_cp_count_base;
    const int64_t *(*count)(const FmiSeed*)          = &fmi_seed_count;
    int64_t        (*sentinel)(const FmiSeed*)       = &fmi_seed_sentinel;
    void           (*sa_prefetch)(FmiSeed*, SMEM*, int64_t*, int64_t*, int64_t,
                                   int32_t, int, int64_t*) = &fmi_seed_sa_prefetch;
    (void)open; (void)close; (void)cp_occ; (void)cp_planes; (void)cp_base; (void)count; (void)sentinel; (void)sa_prefetch;

    // Compile+link smoke only; no real index is opened by a bare invocation.
    // A real prefix argument opts into an actual open/close round-trip.
    if (argc == 2) {
        FmiSeed *h = fmi_seed_open(argv[1]);
        (void)fmi_seed_cp_occ(h);
        (void)fmi_seed_cp_occ_planes(h);
        (void)fmi_seed_cp_count_base(h);
        (void)fmi_seed_count(h);
        (void)fmi_seed_sentinel(h);
        int64_t id = 0;
        fmi_seed_sa_prefetch(h, &s, NULL, NULL, 0, 0, 0, &id);
        fmi_seed_close(h);
    }
    return 0;
}
