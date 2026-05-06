// src/kernel_dispatch.h
//
// Symbol-mangling header for multi-tier kernel builds.
//
// Each kernel TU (`bandedSWA.cpp`, `kswv.cpp`, `ksw.cpp`, `sam_encode.cpp`)
// is compiled multiple times under different `-m...` flags (sse41, sse42,
// avx, avx2, avx512bw on x86_64) with `-DKERNEL_VARIANT=_<suffix>` set on
// the compile line. This header renames every kernel-exported symbol to
// `<name><suffix>` so all per-tier object files can be linked into the same
// binary without symbol collision.
//
// When KERNEL_VARIANT is unset (arm64 build, single-tier x86 fallback),
// names pass through unchanged.
//
// Each kernel TU must include this header *before* any other project
// header that declares the symbols. Kernel headers (bandedSWA.h, kswv.h,
// ksw.h, sam_encode.h) must also include this header at the top so
// declarations and definitions agree.
//
// Non-kernel TUs (bwamem.cpp, fastmap.cpp, etc.) must NEVER define
// KERNEL_VARIANT — they call kernels indirectly through the dispatch
// tables/factories in simd_dispatch.cpp.

#ifndef BWAMEM3_KERNEL_DISPATCH_H
#define BWAMEM3_KERNEL_DISPATCH_H

#ifdef KERNEL_VARIANT
#  define BWAMEM3_PASTE2(a, b) a ## b
#  define BWAMEM3_PASTE(a, b)  BWAMEM3_PASTE2(a, b)

   /* C++ class names — banded SW kernel */
#  define BandedPairWiseSW    BWAMEM3_PASTE(BandedPairWiseSW,    KERNEL_VARIANT)

   /* C++ class names — vectorized KSW batch kernel */
#  define kswv                BWAMEM3_PASTE(kswv,                KERNEL_VARIANT)

   /* C linkage — ksw single-pair entry points (extend, global, align variants).
    * All six are exported from ksw.cpp and called from non-kernel TUs
    * (bwa.cpp, bwamem.cpp, bwamem_pair.cpp). Phase 6 will add C-linkage
    * dispatch wrappers in simd_dispatch.cpp so the unmangled name resolves
    * at link time. */
#  define ksw_extend2         BWAMEM3_PASTE(ksw_extend2,         KERNEL_VARIANT)
#  define ksw_global2         BWAMEM3_PASTE(ksw_global2,         KERNEL_VARIANT)
#  define ksw_extend          BWAMEM3_PASTE(ksw_extend,          KERNEL_VARIANT)
#  define ksw_global          BWAMEM3_PASTE(ksw_global,          KERNEL_VARIANT)
#  define ksw_align2          BWAMEM3_PASTE(ksw_align2,          KERNEL_VARIANT)
#  define ksw_align           BWAMEM3_PASTE(ksw_align,           KERNEL_VARIANT)

   /* C linkage — SAM seq/qual encoder */
#  define sam_encode_seq_fwd  BWAMEM3_PASTE(sam_encode_seq_fwd,  KERNEL_VARIANT)
#  define sam_encode_seq_rev  BWAMEM3_PASTE(sam_encode_seq_rev,  KERNEL_VARIANT)
#  define sam_encode_qual_rev BWAMEM3_PASTE(sam_encode_qual_rev, KERNEL_VARIANT)

   /* C linkage — per-tier kernel factories. Defined in each kernel TU
    * (compiled per-tier on x86); called from simd_dispatch.cpp to construct
    * the right per-tier concrete class without needing the class layout
    * visible at the dispatcher's translation unit. */
#  define make_bsw_kernel    BWAMEM3_PASTE(make_bsw_kernel,    KERNEL_VARIANT)
#  define make_kswv_kernel   BWAMEM3_PASTE(make_kswv_kernel,   KERNEL_VARIANT)
#endif

#endif /* BWAMEM3_KERNEL_DISPATCH_H */
