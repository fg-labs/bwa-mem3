#include "fmi_seed_api.h"
#include <cstdio>
int main() {
    CP_OCC blk; SMEM s; s.k = 0; s.s = 0; (void)blk; (void)s;
    unsigned long m = one_hot_mask_array[3];
    int pc = _mm_countbits_64(m);
    printf("CP_SHIFT=%d CP_MASK=%d popcount=%d\n", CP_SHIFT, CP_MASK, pc);
    return 0;
}
