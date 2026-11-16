// test/framework/x86_tier_caps.h
//
// The CPU check a test must pass before it calls an avx512bw kernel entry point
// directly (make_bsw_kernel_avx512bw, ksw_global2_avx512bw, sam_encode_*_avx512bw,
// ...). Those objects are built with AVX-512 F, BW, VL and DQ (Makefile
// KERNEL_FLAGS_avx512bw), so AVX-512BW alone is not enough: the check asks the
// dispatcher's own kernel-tier ranking, so the tests and production dispatch
// cannot disagree about which hosts run them.
#ifndef BWA_TESTS_X86_TIER_CAPS_H
#define BWA_TESTS_X86_TIER_CAPS_H

#include "simd_dispatch.h"

#if defined(__x86_64__) || defined(__i386__)
namespace bwa_tests {

inline bool host_runs_avx512bw_kernels()
{
    __builtin_cpu_init();
    return bwamem3_x86_kernel_tier(
               __builtin_cpu_supports("avx512bw") != 0, __builtin_cpu_supports("avx512vl") != 0,
               __builtin_cpu_supports("avx512dq") != 0, __builtin_cpu_supports("avx2") != 0,
               __builtin_cpu_supports("avx") != 0, __builtin_cpu_supports("sse4.2") != 0,
               __builtin_cpu_supports("sse4.1") != 0) == BWAMEM3_TIER_AVX512BW;
}

} // namespace bwa_tests
#endif

#endif
