# BWA-MEM3 Apple Silicon Optimization Status

## Current State
- Branch: `apple-silicon`
- Base commit: `7c847ff` - Initial Apple Silicon port complete
- Benchmark baseline (100K high-error reads): ~15.4s (sse2neon), ~14.4s (with kswv NEON)

## Optimization Tasks

| # | Task | Status | Impact | Notes |
|---|------|--------|--------|-------|
| 1 | Correctness verification | ✓ DONE | PASS | 200,006 alignments, 0 differences |
| 2 | Dynamic batch sizing (L2 cache) | ✓ DONE | ~0% | Already at 1024 (compile-time), detection works |
| 3 | Native NEON for bandedSWA.cpp | ✓ DONE | ~4% | Optimized blendv16 with vbsl |
| 4 | Multi-binary launcher | N/A | 0% | Not needed on ARM (single NEON level) |
| 5 | Accelerate.framework usage | ✓ DONE | ~0% | Linked but no suitable compute patterns |
| 6 | M1/M2/M3/M4 detection | ✓ DONE | ~0% | Already in HTStatus(), P/E-cores detected |
| 7 | Native NEON for FMI_search.cpp | N/A | 0% | Memory-bound, not SIMD-friendly |
| 8 | Profile-guided optimization | ✓ DONE | ~3% | Makefile pgo-* targets added |

## Benchmark Results

### Baseline (before optimizations)
- Test: 100K read pairs, 5% error rate, 30% indels, chr17.fa, 8 threads
- sse2neon only: avg 15.4s
- With kswv.cpp NEON: avg 14.4s (~7% faster)

### After Each Optimization

**Task 2 - Dynamic batch sizing:**
- L2 cache detection: 4 MB detected via sysctl
- Dynamic batch size calculation: 1024 (matches compile-time)
- BATCH_SIZE already set to 1024 for Apple Silicon in macro.h
- Platform info now displays at runtime: P-cores, E-cores, L2 cache
- Benchmark: avg 14.9s (similar to baseline, as expected since BATCH_SIZE unchanged)

## Log

### Task 2: Dynamic Batch Sizing
- Added get_l2_cache_size() and get_dynamic_batch_size() to fastmap.cpp
- Fixed HTStatus() not being called on non-NUMA systems (macOS)
- L2 cache correctly detected as 4MB
- Compile-time BATCH_SIZE=1024 is optimal for this cache size
- No performance change expected since batch size was already correct

### Task 3: Native NEON for bandedSWA.cpp
- Optimized _mm_blendv_epi16 to use native NEON vbsl (bitwise select)
- Added `_mm_movemask_epi16` and `_mm_blendv_epi16_fast` to `simd_compat.h`
- Modified bandedSWA.cpp to use vbsl directly on ARM instead of OR/AND/ANDNOT
- Benchmark: 14.4s → 13.8s average (~4% improvement)
- Correctness verified: identical AS:i scores

### Task 6: M1/M2/M3/M4 Detection
- Already implemented in HTStatus() from earlier work
- Detects P-cores and E-cores via hw.perflevel sysctl
- Reports platform info at startup

### Task 8: Profile-Guided Optimization
- Added Makefile targets: pgo-generate, pgo-use, pgo-clean
- Built instrumented binary, ran training workload
- PGO build: avg 14.2s vs baseline 14.6s (~3% improvement)
- Correctness verified: identical AS:i scores
- Modest but measurable improvement

### Task 4: Multi-binary Launcher
- On ARM/Apple Silicon, there's only one NEON instruction set level
- Unlike x86 (SSE41/SSE42/AVX/AVX2/AVX512), no need for multiple binaries
- arm64 Makefile target already creates symlink bwa-mem3 -> bwa-mem3.arm64
- N/A: No benefit from multi-binary approach on ARM

### Task 5: Accelerate.framework
- Framework already linked in build for macOS
- Analyzed compute patterns in bwa-mem3:
  - Smith-Waterman: Already SIMD-optimized (sse2neon)
  - FM-index: Memory-bound, pointer chasing (not vectorizable)
  - Sorting: Small arrays, not suitable for vDSP
- No suitable compute patterns for Accelerate (BLAS/vDSP are for large matrices/vectors)

### Task 7: Native NEON for FMI_search.cpp
- Analyzed backwardExt (85 profile samples) and getSMEMs functions
- Algorithm is fundamentally memory-bound:
  - Sequential dependencies (each step depends on previous result)
  - Pointer chasing through cp_occ arrays
  - Prefetching already implemented
- Popcount (_mm_countbits_64) already uses optimal __builtin_popcountl
- N/A: SIMD won't help memory-bound workload

## Final Summary

### Performance Results
- **Baseline (sse2neon only)**: ~15.4s
- **With kswv.cpp native NEON**: ~14.4s (~7% faster)
- **With bandedSWA.cpp NEON blendv**: ~13.8s (~4% faster)
- **Final optimized build**: ~13.8s average
- **Overall improvement**: ~10% from pure sse2neon baseline

### Key Findings
1. **Native NEON for kswv.cpp**: Most impactful optimization (~7% gain)
2. **Native NEON for bandedSWA.cpp**: Optimized blendv with vbsl (~4% gain)
3. **Memory-bound workloads**: FM-index can't benefit from SIMD
4. **PGO has modest impact**: ~3% improvement, worth using for production builds
5. **Apple Silicon detection**: Works correctly, reports P/E cores and L2 cache

### Recommendations for Production
1. Use `make pgo-generate && <run training> && make pgo-use` for best performance
2. The regular `make arch=arm64` build is adequate for most uses
3. Combined optimizations provide ~10% improvement over baseline sse2neon

### Files Modified
- `src/kswv.cpp`, `src/kswv.h` - Native NEON implementation
- `src/fastmap.cpp` - L2 cache detection, HTStatus for non-NUMA
- `src/macro.h` - BATCH_SIZE tuning for Apple Silicon
- `src/bandedSWA.h` - SIMD width definitions for ARM
- `src/simd_compat.h` - sse2neon integration, memory allocation
- `Makefile` - PGO targets, ARM64 build support

---
*Last updated: All 8 tasks complete*
