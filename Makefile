##/*************************************************************************************
##                           The MIT License
##
##   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
##   Copyright (C) 2019  Intel Corporation, Heng Li.
##
##   Permission is hereby granted, free of charge, to any person obtaining
##   a copy of this software and associated documentation files (the
##   "Software"), to deal in the Software without restriction, including
##   without limitation the rights to use, copy, modify, merge, publish,
##   distribute, sublicense, and/or sell copies of the Software, and to
##   permit persons to whom the Software is furnished to do so, subject to
##   the following conditions:
##
##   The above copyright notice and this permission notice shall be
##   included in all copies or substantial portions of the Software.
##
##   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
##   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
##   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
##   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
##   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
##   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
##   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
##   SOFTWARE.
##
##Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
##                                Heng Li <hli@jimmy.harvard.edu> 
##*****************************************************************************************/

ifneq ($(portable),)
	STATIC_GCC=-static-libgcc -static-libstdc++
endif

EXE=		bwa-mem3
#CXX=		icpc
# Pair CC with CXX for the toolchains we know about so libsais.c (the only
# C TU in this build) doesn't silently fall back to make's default `cc`,
# which can drift from CXX (e.g. CXX=icpx + CC=gcc would mix toolchains).
ifeq ($(CXX), icpc)
	CC= icc
else ifeq ($(CXX), icpx)
	CC= icx
else ifeq ($(CXX), g++)
	CC= gcc
else ifeq ($(CXX), clang++)
	CC= clang
else ifeq ($(CXX), c++)
	# `c++` is GNU make's default and on every system we ship to it is
	# symlinked to the same toolchain as `cc`, so the pairing is safe.
	CC= cc
else
    # Wrappers (`ccache g++`), versioned binaries (`g++-13`), and
    # cross-compilers don't match any case above. Only warn if CC truly
    # is the make-default — operators who pin both CXX and CC explicitly
    # already know what they're doing and shouldn't see a spurious hint.
    ifeq ($(origin CC),default)
        $(warning Unrecognized CXX='$(CXX)'; CC will fall back to default '$(CC)' which may not match CXX. Set CC explicitly to avoid mixing toolchains for libsais.c.)
    endif
endif

# AddressSanitizer support for catching kswv rowMax / SIMD store overruns
# in regression tests (e.g. kswv_nrow_zero_test). Opt-in with `make ASAN=1 ...`
# Forces USE_MIMALLOC off: mimalloc's malloc override interposes before asan
# and the two can't coexist cleanly. Must be set before the mimalloc block so
# USE_MIMALLOC=0 actually takes effect there. CXXFLAGS picks up $(ASAN_FLAGS)
# unconditionally later in the file; it stays empty when ASAN is unset.
ifneq ($(strip $(ASAN)),)
    USE_MIMALLOC = 0
    ASAN_FLAGS   = -fsanitize=address -fno-omit-frame-pointer -O1
    LDFLAGS     += -fsanitize=address
    CFLAGS      += $(ASAN_FLAGS)
endif

# TESTING_BUILD=1 defines BWAMEM3_TESTING which enables test-only injection
# hooks (BWAMEM3_TESTING_HOST_TIER env var read by simd_dispatch.cpp). Used
# by test/regression/host_floor_enforce.sh. Production builds leave this
# unset; integration tests build with `make TESTING_BUILD=1`.
ifneq ($(TESTING_BUILD),)
    CXXFLAGS += -DBWAMEM3_TESTING
endif

# mimalloc integration. Default on — see FG-MAIN.md.
# Override with USE_MIMALLOC=0 to build a stock bwa-mem3 without mimalloc.
USE_MIMALLOC ?= 1

# Detect architecture
UNAME_M := $(shell uname -m)
UNAME_S := $(shell uname -s)
# Treat macOS ("arm64") and Linux ("aarch64") as the same ARM build target.
IS_ARM := $(filter $(UNAME_M),arm64 aarch64)

# Where mimalloc lives and where its CMake build writes artifacts.
MIMALLOC_SRC   = ext/mimalloc
MIMALLOC_BUILD = $(MIMALLOC_SRC)/build

# Per-platform library basename. Linux: static archive. macOS: dynamic lib
# (mimalloc's malloc override on macOS requires a dylib + dyld interposing).
ifeq ($(UNAME_S),Darwin)
    MIMALLOC_LIB = $(MIMALLOC_BUILD)/libmimalloc.dylib
    MIMALLOC_CMAKE_FLAGS = -DMI_BUILD_SHARED=ON -DMI_BUILD_STATIC=OFF \
                           -DMI_BUILD_OBJECT=OFF -DMI_BUILD_TESTS=OFF \
                           -DMI_OVERRIDE=ON -DCMAKE_BUILD_TYPE=Release
else
    MIMALLOC_LIB = $(MIMALLOC_BUILD)/libmimalloc.a
    MIMALLOC_CMAKE_FLAGS = -DMI_BUILD_SHARED=OFF -DMI_BUILD_STATIC=ON \
                           -DMI_BUILD_OBJECT=OFF -DMI_BUILD_TESTS=OFF \
                           -DMI_OVERRIDE=ON -DCMAKE_BUILD_TYPE=Release
endif

# Link flags that inject mimalloc's malloc overrides.
# On Linux, --whole-archive forces the linker to keep symbols it would
# otherwise drop (malloc/free would come from libc first). On macOS, we
# just link the dylib; dyld interposes at load time — no --whole-archive
# equivalent is needed (and -force_load does NOT enable malloc interpose).
#
# On macOS we set two rpaths: @executable_path/. is the portable one used
# when the binary ships alongside libmimalloc.dylib; the $(abspath ...)
# rpath is a dev-only fallback that lets the binary run in-tree without
# first copying the dylib. For distribution, the @executable_path rpath
# resolves first and the abspath is harmlessly ignored (or can be removed
# with `install_name_tool -delete_rpath`).
ifeq ($(USE_MIMALLOC),1)
    ifeq ($(UNAME_S),Darwin)
        MIMALLOC_LDFLAGS = -L$(MIMALLOC_BUILD) -lmimalloc \
                           -Wl,-rpath,@executable_path/. \
                           -Wl,-rpath,$(abspath $(MIMALLOC_BUILD))
    else
        MIMALLOC_LDFLAGS = -Wl,--whole-archive $(MIMALLOC_LIB) -Wl,--no-whole-archive
    endif
    CPPFLAGS += -DUSE_MIMALLOC=1
else
    MIMALLOC_LDFLAGS =
endif

# ARM/Apple Silicon support
ifneq ($(IS_ARM),)
    ARCH_FLAGS = -DAPPLE_SILICON=1
    # sse2neon flags - define SSE feature macros for translation
    SSE2NEON_FLAGS = -D__SSE__=1 -D__SSE2__=1 -D__SSE3__=1 -D__SSSE3__=1 -D__SSE4_1__=1 -D__SSE4_2__=1
    SSE2NEON_INCLUDES = -Iext/sse2neon
    CPPFLAGS += $(SSE2NEON_FLAGS)
    INCLUDES += $(SSE2NEON_INCLUDES)
    # Apple Silicon uses 128-byte cache lines
    CPPFLAGS += -DCACHE_LINE_BYTES=128
    # Link Accelerate framework on macOS for potential BLAS/vecLib usage
    ifeq ($(UNAME_S),Darwin)
        LIBS_EXTRA = -framework Accelerate
    endif
else
    ARCH_FLAGS = -msse -msse2 -msse3 -mssse3 -msse4.1
endif

CPPFLAGS+=	-DENABLE_PREFETCH -DV17=1 -DMATE_SORT=0 -DLIBSAIS_OPENMP

# Version string for `bwa-mem3 version` and the @PG VN: field. Prefer
# `git describe` (e.g. v2.3-30-g61813ef, with -dirty suffix for modified
# trees) so the stamped version always reflects the actual build. Fall
# back to a static tag for source-tarball / shallow-clone builds.
FG_LABS_VERSION_FALLBACK := 0.2.0
VERSION_STRING := $(shell git describe --tags --dirty 2>/dev/null || echo $(FG_LABS_VERSION_FALLBACK))
INCLUDES+=   -Isrc -Iext/safestringlib/include -Iext/htslib -Iext/libsais/include
ifeq ($(USE_MIMALLOC),1)
    INCLUDES += -Iext/mimalloc/include
endif

# libsais (pinned in ext/libsais; see submodule SHA): linear-time suffix
# array / BWT construction via SA-IS. Compiled with OpenMP so
# libsais_gsa_omp can run parallel induced-sorting. libomp is already a
# link dep of bwa-mem3's alignment paths; no new dep.
LIBSAIS_DIR    = ext/libsais
LIBSAIS_OBJS   = $(LIBSAIS_DIR)/src/libsais.o $(LIBSAIS_DIR)/src/libsais64.o
LIBSAIS_CFLAGS = -O3 -std=c99 -DLIBSAIS_OPENMP -I$(LIBSAIS_DIR)/include
ifeq ($(UNAME_S),Darwin)
    # Resolve at parse time but defer the missing-libomp error until the
    # libsais recipe actually runs, so `make clean`, `make print-mimalloc-config`,
    # etc. still work on hosts without libomp installed. The check below in
    # the libsais pattern rule produces the actionable hint when needed.
    LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
    LIBSAIS_OPENMP_CFLAGS = -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
    LIBSAIS_OPENMP_LIBS   = -L$(LIBOMP_PREFIX)/lib -lomp
else
    LIBSAIS_OPENMP_CFLAGS = -fopenmp
    LIBSAIS_OPENMP_LIBS   = -fopenmp
endif

LIBS=		-lpthread -lm -lz -L. -lbwa -Lext/safestringlib -lsafestring -Lext/htslib -lhts $(LIBSAIS_OPENMP_LIBS) $(STATIC_GCC) $(LIBS_EXTRA)
# Non-kernel objects: always compiled once at the baseline ISA and linked into
# libbwa.a on every build (arm64 and x86 alike).
OBJS=		src/fastmap.o src/bwtindex.o src/utils.o src/memcpy_bwamem.o src/kthread.o \
			src/kstring.o src/bntseq.o src/bwamem.o src/profiling.o \
			src/FMI_search.o src/read_index_ele.o src/bwamem_pair.o src/bwa.o \
			src/bwamem_extra.o src/kopen.o src/bam_writer.o src/meth_bam.o \
			src/meth_orig_ref.o src/meth_xm.o \
			src/packed_text.o src/fm_index_writer.o src/index_prelude.o \
			src/system.o src/libsais_build.o \
			src/bwa_shm.o src/simd_dispatch.o

# Kernel TUs (bandedSWA, kswv, ksw, sam_encode) are compiled per-tier on x86
# and linked directly via KERNEL_TIER_OBJS_LINK. The dispatch wrappers in
# simd_dispatch.cpp provide the unmangled entry points. On arm64 there is only
# one tier (NEON, unmangled), so the baseline objects ARE the only copies.
# On x86, the baseline (unsuffixed) objects provide the concrete kswv/ksw/etc.
# class bodies needed by test binaries (e.g. kswv_nrow_zero_test) that
# instantiate these classes directly rather than going through the dispatcher.
KERNEL_BASELINE_OBJS = src/bandedSWA.o src/kswv.o src/ksw.o src/sam_encode.o
OBJS += $(KERNEL_BASELINE_OBJS)
BWA_LIB=    libbwa.a
SAFE_STR_LIB=    ext/safestringlib/libsafestring.a
HTS_LIB=    ext/htslib/libhts.a

# Architecture-specific builds (x86 only, ARM uses default from above)
ifeq ($(IS_ARM),)
ifeq ($(arch),sse41)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-msse4.1
	else
		ARCH_FLAGS=-msse -msse2 -msse3 -mssse3 -msse4.1
	endif
else ifeq ($(arch),sse42)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-msse4.2
	else
		ARCH_FLAGS=-msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2
	endif
else ifeq ($(arch),avx)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-mavx ##-xAVX
	else
		ARCH_FLAGS=-mavx
	endif
else ifeq ($(arch),avx2)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-march=core-avx2 #-xCORE-AVX2
	else
		ARCH_FLAGS=-mavx2
	endif
else ifeq ($(arch),avx512)
	# Legacy alias for arch=avx512bw (preserved for backward compat with
	# pre-PR #16 invocations). Keep flags identical to the avx512bw branch
	# below — including the -mprefer-vector-width=256 / -qopt-zmm-usage=low
	# autovec cap. See the avx512bw branch for the rationale.
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-xCORE-AVX512 -qopt-zmm-usage=low
	else
		ARCH_FLAGS=-mavx512f -mavx512bw -mprefer-vector-width=256
	endif
else ifeq ($(arch),avx512bw)
	# Explicit BW target: double the lane width vs AVX2 (64x8-bit / 32x16-bit).
	# AVX-512BW implies AVX-512F; -mavx512bw alone enables BW+F on gcc/clang
	# but we list both flags for clarity.
	#
	# -mprefer-vector-width=256 (gcc/clang) / -qopt-zmm-usage=low (icpc):
	# keep AVX-512BW *capabilities* available (32 zmm registers, mask
	# registers, byte/word lane permutes, gather/scatter) but cap the
	# auto-vectorizer's preferred SIMD width at 256-bit. Two effects:
	#   1. AMD Zen 4 (c7a / Genoa) splits 512-bit AVX-512 ops into
	#      2x 256-bit µops per op. For short-trip-count auto-vec loops
	#      that's a regression — 2x latency without amortizing the
	#      reduced iteration count, plus more I-cache pressure. Capping
	#      at 256-bit keeps the compiler from widening those loops.
	#   2. Intel Sapphire Rapids (c7i / m7i) has native 512-bit
	#      execution but pays a ~3-5% AVX-512 frequency downclock under
	#      sustained heavy use, plus AVX-512↔AVX2 transition penalties
	#      when non-kernel TUs running 512-bit code call into
	#      explicitly-256-bit kernel TUs. Capping at 256-bit avoids
	#      both.
	# Empirical: c7a wgs-5M shm-warmed -4.4% wall vs avx2 baseline
	# (vanilla avx512bw was -2.2%, the cap adds another 2% on top).
	# c7i wgs-5M is a wash either way (-0.7%). The hand-tuned 512-bit
	# kernel TUs are unaffected — they use intrinsics, not auto-vec.
	# Canonical mitigation per FFmpeg, libvpx, ISPC docs.
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-xCORE-AVX512 -qopt-zmm-usage=low
	else
		ARCH_FLAGS=-mavx512f -mavx512bw -mprefer-vector-width=256
	endif
else ifeq ($(arch),native)
	ARCH_FLAGS=-march=native
else ifneq ($(arch),)
# To provide a different architecture flag like -march=core-avx2.
	ARCH_FLAGS=$(arch)
else
myall:single
endif
endif

# ARM64/Apple Silicon single-binary build
ifneq ($(IS_ARM),)
ifeq ($(arch),arm64)
    ARCH_FLAGS = -DAPPLE_SILICON=1
else ifeq ($(arch),)
myall:arm64
endif
endif

CXXFLAGS+=	-g -O3 -std=gnu++14 -fpermissive $(ARCH_FLAGS) $(ASAN_FLAGS) $(LIBSAIS_OPENMP_CFLAGS) $(EXTRA_CXXFLAGS) #-Wall ##-xSSE2

# CXXFLAGS used for per-tier kernel TU compilation. Must NOT contain any
# -m... ISA flag — each per-tier rule appends the right one.
BASE_CXXFLAGS = -g -O3 -std=gnu++14 -fpermissive $(ASAN_FLAGS) $(LIBSAIS_OPENMP_CFLAGS) $(EXTRA_CXXFLAGS)

# Per-tier ISA flag groups for kernel multi-tier compilation (x86_64 only).
KERNEL_FLAGS_sse41    = -msse4.1
KERNEL_FLAGS_sse42    = -msse4.2 -mpopcnt
KERNEL_FLAGS_avx      = -mavx -mpopcnt
KERNEL_FLAGS_avx2     = -mavx -mavx2 -mpopcnt
KERNEL_FLAGS_avx512bw = -mavx -mavx2 -mavx512f -mavx512bw -mpopcnt

# Source files compiled at every tier on x86_64.
KERNEL_SRCS = src/bandedSWA.cpp src/kswv.cpp src/ksw.cpp src/sam_encode.cpp

# All per-tier kernel objects. On x86_64 these are 5 tiers × 4 kernel TUs = 20 .o files.
ifneq ($(IS_ARM),)
    # arm64: single unsuffixed build (KERNEL_VARIANT unset). Already in $(OBJS).
    KERNEL_TIER_OBJS =
else
    KERNEL_TIER_OBJS = \
        $(patsubst src/%.cpp,src/%.sse41.o,$(KERNEL_SRCS)) \
        $(patsubst src/%.cpp,src/%.sse42.o,$(KERNEL_SRCS)) \
        $(patsubst src/%.cpp,src/%.avx.o,$(KERNEL_SRCS)) \
        $(patsubst src/%.cpp,src/%.avx2.o,$(KERNEL_SRCS)) \
        $(patsubst src/%.cpp,src/%.avx512bw.o,$(KERNEL_SRCS))
endif

# COVERAGE=1 augments CXXFLAGS/LDFLAGS with --coverage and overrides -O3 with
# -O0 so gcov line numbers correspond 1:1 with source. Consumed by the CI
# `coverage` job; not part of any shipped binary.
ifneq ($(COVERAGE),)
    CXXFLAGS += -O0 --coverage
    LDFLAGS  += --coverage
endif

# Control build flag for the batched mate-rescue SW port on ARM.
# When set (e.g. `make arm64 DISABLE_BATCHED_MATESW=1`), the source gate for
# the new batched path falls through to the legacy scalar mem_sam_pe. Used by
# the proto-neon-kswv CI to A/B the same commit with the port on vs. off.
# Pass the caller-supplied value through verbatim so `DISABLE_BATCHED_MATESW=0`
# still selects the batched path (ifdef would be true even for =0).
ifneq ($(strip $(DISABLE_BATCHED_MATESW)),)
    CPPFLAGS += -DDISABLE_BATCHED_MATESW=$(DISABLE_BATCHED_MATESW)
endif

.PHONY:all myall arm64 clean depend single all-single print-mimalloc-config kswv_nrow_zero_test shm_section_find_test shm_pack_round_trip_test shm_lock_destroy_test test test-injection FORCE pgo-generate pgo-use pgo-clean profile-build profile-clean lto-build lto-clean docs docs-serve docs-cli docs-clean docs-install-tools
.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

# Per-tier kernel compile rules. Active only on x86_64 multi-tier builds.
src/%.sse41.o: src/%.cpp
	$(CXX) -c $(BASE_CXXFLAGS) $(KERNEL_FLAGS_sse41) $(CPPFLAGS) -DKERNEL_VARIANT=_sse41 $(INCLUDES) $< -o $@

src/%.sse42.o: src/%.cpp
	$(CXX) -c $(BASE_CXXFLAGS) $(KERNEL_FLAGS_sse42) $(CPPFLAGS) -DKERNEL_VARIANT=_sse42 $(INCLUDES) $< -o $@

src/%.avx.o: src/%.cpp
	$(CXX) -c $(BASE_CXXFLAGS) $(KERNEL_FLAGS_avx) $(CPPFLAGS) -DKERNEL_VARIANT=_avx $(INCLUDES) $< -o $@

src/%.avx2.o: src/%.cpp
	$(CXX) -c $(BASE_CXXFLAGS) $(KERNEL_FLAGS_avx2) $(CPPFLAGS) -DKERNEL_VARIANT=_avx2 $(INCLUDES) $< -o $@

src/%.avx512bw.o: src/%.cpp
	$(CXX) -c $(BASE_CXXFLAGS) $(KERNEL_FLAGS_avx512bw) $(CPPFLAGS) -DKERNEL_VARIANT=_avx512bw $(INCLUDES) $< -o $@

all:$(EXE)

# Note: simd_dispatch.cpp is compiled at the regular BASELINE_ARCH like
# every other non-kernel TU. An earlier draft of this PR compiled it at
# -march=x86-64 to keep the precheck path SIGILL-safe on too-old hosts,
# but that broke `g_build_tier` (a `static constexpr` in this TU derived
# from __AVX2__/__SSE4_1__/etc., which only get defined when the matching
# -m flag is in scope). With the override, every binary reported its
# floor as "scalar" and the precheck silently became a no-op.
#
# In practice the precheck path is scalar-only — std::call_once,
# integer comparisons, getenv, snprintf, fputs, exit — with no array
# loops the compiler could autovectorize. main.cpp's preamble before
# the precheck (argc check, rdtsc, sleep, argv strcmp scan) is the
# same shape: scalar branching + libc calls. Even at -mavx2 the
# compiler doesn't emit AVX2 for any of this, so the precheck actually
# fires on EVERY host below the build floor — sse42, sse41, scalar,
# cross-family — not just one tier below. Verified empirically with
# the BWAMEM3_TESTING_HOST_TIER injection test.

# Regenerate src/version.h on every invocation, but only touch the file
# (and thus trigger a main.o rebuild) when the string actually changed.
# Must be declared after `all:$(EXE)` so FORCE is never picked as the
# default goal when the caller supplies `arch=...` (which skips the
# `myall:` dispatch branch).
FORCE:
src/version.h: FORCE
	@printf '#ifndef BWA_MEM3_VERSION_H\n#define BWA_MEM3_VERSION_H\n#define PACKAGE_VERSION "%s"\n#endif\n' '$(VERSION_STRING)' > $@.tmp
	@if ! cmp -s $@.tmp $@ 2>/dev/null; then mv $@.tmp $@; else rm -f $@.tmp; fi

src/main.o: src/version.h

# Baseline ISA tier for non-kernel TUs in the x86 single-binary build.
# Defaults to avx2: every host that runs bwa-mem3 in practice has AVX2
# (Haswell, 2013+; any host with AVX-512 also has AVX2), and dropping the
# baseline below avx2 measurably slows hot non-kernel paths (chain
# extension, FMI BWT walks, mate scoring) because the compiler can no
# longer auto-vectorize them at 256-bit width. Override to sse41 (or
# sse42, avx) for vintage hardware; the per-tier kernel objects are still
# compiled at every tier regardless, so kernel dispatch on lower-tier
# hosts continues to work.
BASELINE_ARCH ?= avx2

# Single-binary multi-tier build. All kernel TUs are compiled at every
# tier; the dispatcher picks the right per-tier subclass at runtime.
# Replaces the `multi` target's 5 sequential clean rebuilds + execv launcher.
.PHONY: single
single: $(if $(filter 1,$(USE_MIMALLOC)),$(MIMALLOC_LIB))
ifneq ($(IS_ARM),)
	@echo "ARM64 detected - building single arm64 binary instead of multi-tier"
	$(MAKE) arm64
else
	$(MAKE) arch=$(BASELINE_ARCH) EXE=bwa-mem3 CXX="$(CXX)" KERNEL_TIER_OBJS_LINK="$(KERNEL_TIER_OBJS)" all-single
endif

# Internal: builds the single binary with the BASELINE_ARCH tier for
# non-kernel TUs and links all KERNEL_TIER_OBJS_LINK on top.
.PHONY: all-single
all-single: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) $(LIBSAIS_OBJS) $(KERNEL_TIER_OBJS_LINK) src/main.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) src/main.o $(KERNEL_TIER_OBJS_LINK) $(BWA_LIB) $(LIBSAIS_OBJS) $(LIBS) $(MIMALLOC_LDFLAGS) -o $(EXE)

# ARM64/Apple Silicon build target - single binary, no multi-binary launcher needed
arm64:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=arm64 EXE=bwa-mem3.arm64 CXX="$(CXX)" all
	ln -sf bwa-mem3.arm64 bwa-mem3


$(EXE):$(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) $(LIBSAIS_OBJS) $(if $(filter 1,$(USE_MIMALLOC)),$(MIMALLOC_LIB)) src/main.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) src/main.o $(BWA_LIB) $(LIBSAIS_OBJS) $(LIBS) $(MIMALLOC_LDFLAGS) -o $@

# Regression test for issue 38 / upstream PR 289: exercises an all-len1==0
# batch that drives each SIMD kswv kernel through the nrow==0 path. Without
# the post-loop `if (i > 0)` guard, the rowMax store writes SIMD_WIDTH* bytes
# before the allocation and aborts at a later allocator operation; under
# asan the write is reported directly.
#
# On x86 multi-tier builds, libbwa.a's baseline kswv.o is compiled at the
# BASELINE_ARCH tier (avx2 by default; sse41 if overridden, in which case
# the SSE-only stub that calls exit() would fire). Compile a separate
# native-tier copy of kswv.cpp (src/kswv.native.o) and link it ahead of
# libbwa.a so the linker picks the host's native-ISA concrete kswv class
# regardless of BASELINE_ARCH. On arm64 -march=native resolves to the NEON
# path already covered by the baseline objects.
src/kswv.native.o: src/kswv.cpp
	$(CXX) -c $(BASE_CXXFLAGS) -march=native $(CPPFLAGS) $(INCLUDES) $< -o $@

kswv_nrow_zero_test: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) src/kswv.native.o test/kswv_nrow_zero_test.o
	$(CXX) $(BASE_CXXFLAGS) -march=native $(LDFLAGS) test/kswv_nrow_zero_test.o src/kswv.native.o $(BWA_LIB) $(LIBS) -o $@

# Build the test binaries with the same ARCH_FLAGS as libbwa.a so the
# test binary's kswv.h preprocessor state (SIMD_WIDTH8, BWA_TESTS_HAVE_KSWV)
# matches what libbwa.a was compiled with. Consumed by ci.yml so that e.g.
# the sse41 matrix row builds test/framework with -msse4.1 only (matching
# libbwa.a, which then lacks kswv::getScores8 — the BWA_TESTS_HAVE_KSWV
# macro guards the test away).
.PHONY: test-binaries
# $(SAFE_STR_LIB) and $(HTS_LIB) are real link-time deps: test/Makefile's
# bwa_mem3_tests_unit recipe references ../ext/safestringlib/libsafestring.a
# and ../ext/htslib/libhts.a directly. Without these prereqs, callers that
# skip the bwa-mem3 binary build (which builds them as a side-effect of
# $(EXE) deps) link-fail.
test-binaries: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB)
	$(MAKE) -C test framework unit integration \
	    CXX="$(CXX)" \
	    COVERAGE=$(COVERAGE) \
	    ARCH_FLAGS_FROM_PARENT='$(ARCH_FLAGS)'

shm_section_find_test: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) $(LIBSAIS_OBJS) test/shm_section_find_test.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) test/shm_section_find_test.o $(BWA_LIB) $(LIBSAIS_OBJS) $(LIBS) -o $@

shm_pack_round_trip_test: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) $(LIBSAIS_OBJS) test/shm_pack_round_trip_test.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) test/shm_pack_round_trip_test.o $(BWA_LIB) $(LIBSAIS_OBJS) $(LIBS) -o $@

shm_lock_destroy_test: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) $(LIBSAIS_OBJS) test/shm_lock_destroy_test.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) test/shm_lock_destroy_test.o $(BWA_LIB) $(LIBSAIS_OBJS) $(LIBS) -o $@

test/shm_pack_round_trip_test.o: test/shm_pack_round_trip_test.cpp

# Run the in-tree tests via the unit-test harness in test/, plus the
# standalone regressions (kswv_nrow_zero_test + shm_section_find_test +
# shm_lock_destroy_test). shm_pack_round_trip_test runs via
# test/shm_pack_round_trip_test.sh which builds the phiX index first;
# invoked from test/run_unit_tests.sh.
#
# Note: depends on `bwa-mem3` so version_banner.sh has a binary to grep —
# previously `test:` only built the test harness binaries, not the main
# executable.
test: test-binaries kswv_nrow_zero_test shm_section_find_test shm_lock_destroy_test bwa-mem3
	./test/bwa_mem3_tests_unit
	./test/bwa_mem3_tests_integration
	./kswv_nrow_zero_test
	./shm_section_find_test
	./shm_lock_destroy_test
	BWA_MEM3=./bwa-mem3 ./test/regression/version_banner.sh

# Regression test that requires a binary built with TESTING_BUILD=1
# (enables BWAMEM3_TESTING_HOST_TIER env-var injection). Not invoked by
# the default `test` target because it requires a non-production build.
# CI runs `make clean && make TESTING_BUILD=1 && make test-injection`
# as a separate matrix row.
test-injection: bwa-mem3
	BWA_MEM3_TESTING=./bwa-mem3 INJECTED_TIER=sse41 PARITY_FA=/dev/null \
		./test/regression/host_floor_enforce.sh

test/kswv_nrow_zero_test.o: test/kswv_nrow_zero_test.cpp
	$(CXX) -c $(BASE_CXXFLAGS) -march=native $(CPPFLAGS) $(INCLUDES) $< -o $@

test/shm_section_find_test.o: test/shm_section_find_test.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

test/shm_lock_destroy_test.o: test/shm_lock_destroy_test.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

# Archive both the baseline (unmangled) kernel objects from $(OBJS) and the
# per-tier (mangled) objects from $(KERNEL_TIER_OBJS) into libbwa.a. The
# duplication with $(KERNEL_TIER_OBJS_LINK) on the all-single recipe is
# intentional: legacy single-arch builds (`make arch=avx2`, `make arch=sse41`,
# etc.) link only libbwa.a on the binary line and rely on libbwa.a to satisfy
# simd_dispatch.o's references to make_kswv_kernel_<tier>, ksw_*_<tier>, and
# sam_encode_*_<tier>. Dropping per-tier objects from libbwa.a breaks those
# matrix rows. On arm64 $(KERNEL_TIER_OBJS) is empty.
$(BWA_LIB):$(OBJS) $(KERNEL_TIER_OBJS)
	ar rcs $(BWA_LIB) $(OBJS) $(KERNEL_TIER_OBJS)

# safestringlib's safeclib_private.h calls abort() / memcpy() via macros
# without including <stdlib.h> / <string.h> in the TU, which fails the
# implicit-function-declaration check in clang >= 15. Some safeclib TUs
# (strcasecmp_s.c, strcasestr_s.c) likewise call toupper() without
# including <ctype.h>. Force-include stdlib.h + ctype.h whenever the
# build CC is clang (Linux clang job) or whenever we're on Darwin (where
# the system `cc` is always clang). On Darwin we also redefine memset_s:
# macOS libc declares C11 Annex K memset_s with a different signature,
# which conflicts with the safestringlib definition.
SAFE_EXTRA_CFLAGS =
SAFE_CC_BASENAME := $(notdir $(CC))
ifeq ($(UNAME_S),Darwin)
    SAFE_EXTRA_CFLAGS += -include stdlib.h -include ctype.h -Dmemset_s=_safestringlib_memset_s
else ifneq (,$(findstring clang,$(SAFE_CC_BASENAME)))
    SAFE_EXTRA_CFLAGS += -include stdlib.h -include ctype.h
endif

$(SAFE_STR_LIB):
	cd ext/safestringlib/ && $(MAKE) clean && $(MAKE) CC="$(CC)" CFLAGS="-Iinclude -Isafeclib $(SAFE_EXTRA_CFLAGS) -fstack-protector-strong -fPIE -fPIC -O2" directories libsafestring.a

# htslib: minimal configure (no lzma/bz2/curl/S3/GCS/plugins), zlib only.
# Guard on config.mk (only created by ./configure) rather than Makefile, which
# is checked into the htslib tree and would make the guard a no-op.
$(HTS_LIB):
	cd ext/htslib && \
	    ([ -f config.mk ] || (autoreconf -i && \
	        ./configure --disable-lzma --disable-libcurl --disable-gcs \
	                    --disable-s3 --disable-plugins --disable-bz2)) && \
	    $(MAKE) libhts.a

# libsais: compile the two C sources we use (libsais.c + libsais64.c) as
# plain .o files. OpenMP enabled via LIBSAIS_OPENMP so libsais64_gsa_omp
# can run parallel induced-sorting.
#
# CFLAGS / CPPFLAGS / ASAN_FLAGS are forwarded so that ASAN builds
# instrument libsais and so that package-manager-supplied flags
# (Conda/Homebrew/distro) reach the libsais TU. LIBSAIS_CFLAGS is appended
# last so its -O3/-std=c99 take precedence over any user override.
$(LIBSAIS_DIR)/src/%.o: $(LIBSAIS_DIR)/src/%.c
	@if [ ! -f $(LIBSAIS_DIR)/include/libsais64.h ]; then \
	    echo "ERROR: $(LIBSAIS_DIR) is empty. Run: git submodule update --init --recursive"; \
	    exit 1; \
	fi
	@if [ "$(UNAME_S)" = "Darwin" ] && [ -z "$(strip $(LIBOMP_PREFIX))" ]; then \
	    echo "ERROR: libomp not found. Install with 'brew install libomp', or set LIBOMP_PREFIX to its install prefix."; \
	    exit 1; \
	fi
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $(ASAN_FLAGS) $(LIBSAIS_CFLAGS) $(LIBSAIS_OPENMP_CFLAGS) $< -o $@

# Build mimalloc via its own CMake system. Shells out to cmake once and
# caches the build tree under ext/mimalloc/build. This rule always builds
# when invoked; USE_MIMALLOC=0 consumers simply don't depend on it.
$(MIMALLOC_LIB):
	@if [ ! -f $(MIMALLOC_SRC)/CMakeLists.txt ]; then \
		echo "ERROR: $(MIMALLOC_SRC) is empty. Run: git submodule update --init --recursive"; \
		exit 1; \
	fi
	mkdir -p $(MIMALLOC_BUILD)
	cd $(MIMALLOC_BUILD) && cmake $(MIMALLOC_CMAKE_FLAGS) .. && $(MAKE)

clean: pgo-clean profile-clean lto-clean
	rm -fr src/*.o src/version.h test/*.o $(BWA_LIB) $(EXE) kswv_nrow_zero_test shm_section_find_test shm_pack_round_trip_test shm_lock_destroy_test bwa-mem3.arm64
	rm -f $(LIBSAIS_OBJS)
	rm -f src/*.gcno src/*.gcda
	$(MAKE) -C test clean
	cd ext/safestringlib/ && $(MAKE) clean
	-[ -f ext/htslib/config.mk ] && cd ext/htslib && $(MAKE) distclean
	rm -rf $(MIMALLOC_BUILD)

# ----------------------------------------------------------------------------
# Documentation (mdbook). See docs/superpowers/specs/ for design.
# ----------------------------------------------------------------------------

# Sub-commands whose --help is captured into docs/_generated/cli/.
DOCS_CLI_SUBCMDS := index mem shm version

docs:
	cd docs && mdbook build

docs-serve:
	cd docs && mdbook serve --open

docs-cli: $(EXE)
	@mkdir -p docs/_generated/cli
	@for sub in $(DOCS_CLI_SUBCMDS); do \
		echo "  CAPTURE  docs/_generated/cli/$$sub.txt"; \
		./$(EXE) $$sub --help 2>&1 \
			| sed -e 's/[[:space:]]*$$//' \
			| grep -v '^Total time taken:' \
			| grep -v '^Looking to launch ' \
			| grep -v '^Launching executable ' \
			| grep -v '^SIMD floor: ' \
			| grep -v '^SIMD runtime: ' \
			| grep -v '^\[W::' \
			| awk '/^v?[0-9]+\.[0-9]+/ {print "v<MAJOR.MINOR>-<N>-g<COMMIT>"; next} {print}' \
			> docs/_generated/cli/$$sub.txt; \
	done

docs-clean:
	rm -rf docs/book

docs-install-tools:
	cargo install mdbook --version 0.5.2 --locked
	cargo install mdbook-mermaid --version 0.17.0 --locked

# Profile-Guided Optimization (PGO) targets.
#
# Usage (host-default arch, single shared profile dir — preserves prior
# arm64 behavior on Apple Silicon / aarch64 hosts):
#   make pgo-generate && <run training workload> && make pgo-use
#
# Multi-arch / multi-regime usage (override at command line):
#   make pgo-generate PGO_ARCH=avx2 PGO_PROFILE_DIR=/path/to/regimeA
#   <run training>
#   make pgo-use PGO_ARCH=avx2 PGO_PROFILE_DIR=/path/to/regimeA
#
# PGO_ARCH accepts the same values as the top-level `arch=` knob: arm64,
# sse41, sse42, avx, avx2, avx512, avx512bw, native, or any custom flag
# string. Defaults match the host: arm64 on Apple Silicon / aarch64,
# native otherwise. Output binaries are arch-suffixed when PGO_ARCH is
# non-default, so multiple per-arch builds coexist:
#   PGO_ARCH=arm64  -> bwa-mem3.pgo-instr,    bwa-mem3.pgo
#   PGO_ARCH=avx2   -> bwa-mem3.pgo-instr.avx2, bwa-mem3.pgo.avx2
ifneq ($(IS_ARM),)
    PGO_ARCH ?= arm64
else
    PGO_ARCH ?= native
endif
PGO_PROFILE_DIR ?= pgo_profiles

# Output names: keep the bare names when PGO_ARCH is the default arm64
# (backward-compat); arch-suffix otherwise so per-arch outputs don't collide.
ifeq ($(PGO_ARCH),arm64)
    PGO_INSTR_EXE = bwa-mem3.pgo-instr
    PGO_FINAL_EXE = bwa-mem3.pgo
else
    PGO_INSTR_EXE = bwa-mem3.pgo-instr.$(PGO_ARCH)
    PGO_FINAL_EXE = bwa-mem3.pgo.$(PGO_ARCH)
endif

pgo-generate:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=$(PGO_ARCH) EXE=$(PGO_INSTR_EXE) EXTRA_CXXFLAGS="-fprofile-generate=$(PGO_PROFILE_DIR)" CXX="$(CXX)" all
	@echo "PGO instrumented binary built: $(PGO_INSTR_EXE) (arch=$(PGO_ARCH), profile dir=$(PGO_PROFILE_DIR))"
	@echo "Run training workload with $(PGO_INSTR_EXE), then: make pgo-use PGO_ARCH=$(PGO_ARCH) PGO_PROFILE_DIR=$(PGO_PROFILE_DIR)"

pgo-use:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=$(PGO_ARCH) EXE=$(PGO_FINAL_EXE) EXTRA_CXXFLAGS="-fprofile-use=$(PGO_PROFILE_DIR) -fprofile-correction" CXX="$(CXX)" all
	@echo "PGO optimized binary built: $(PGO_FINAL_EXE) (arch=$(PGO_ARCH))"

pgo-clean:
	rm -rf $(PGO_PROFILE_DIR) bwa-mem3.pgo-instr bwa-mem3.pgo bwa-mem3.pgo-instr.* bwa-mem3.pgo.*

# profile-build / lto-build target arch. Mirrors PGO_ARCH: defaults to
# arm64 on Apple Silicon / aarch64 hosts (preserves prior behavior), and
# native otherwise. Override at the command line for cross-builds, e.g.
#   make profile-build PROFILE_ARCH=avx2
#   make lto-build LTO_ARCH=avx512bw
ifneq ($(IS_ARM),)
    PROFILE_ARCH ?= arm64
    LTO_ARCH     ?= arm64
else
    PROFILE_ARCH ?= native
    LTO_ARCH     ?= native
endif

# Compute-only profile build. -DDISABLE_OUTPUT short-circuits BAM/SAM
# per-record writes AND writer open + header emit, so wall-clock measurements
# exclude all output I/O (no -o file open, no @HD/@SQ/@PG emission). All
# upstream alignment work runs unchanged; the per-stage tprof[] counters
# (printed at end of run) are unaffected.
# Usage: make profile-build
#        make profile-build PROFILE_ARCH=avx2     # cross-build
#        ./bwa-mem3.profile mem -t N idx r1.fq.gz r2.fq.gz
profile-build:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=$(PROFILE_ARCH) EXE=bwa-mem3.profile EXTRA_CXXFLAGS="$(EXTRA_CXXFLAGS) -DDISABLE_OUTPUT" CXX="$(CXX)" all
	@echo "Compute-only profile binary: bwa-mem3.profile (arch=$(PROFILE_ARCH), output I/O skipped)"
	# Drop variant-flagged objects from the shared cache so a subsequent
	# `make all` doesn't relink stale -DDISABLE_OUTPUT objects.
	rm -f src/*.o $(BWA_LIB)

profile-clean:
	rm -f bwa-mem3.profile

# Link-Time Optimization build.
# Usage: make lto-build
#        make lto-build LTO_ARCH=avx2             # cross-build
#        ./bwa-mem3.lto mem -t N idx r1.fq.gz r2.fq.gz
# Compiles all bwa-mem3 sources with LTO and links with LTO. Non-bwa-mem3
# deps (htslib, mimalloc, safestringlib) keep their non-LTO objects; the
# linker still does LTO across bwa-mem3's own .o. On GCC,
# -fno-semantic-interposition additionally allows more aggressive inlining
# across translation units (no effect on clang, silently ignored).

# LTO_FLAG is detected at recipe-time (not Makefile-parse time) so a stale
# or missing $(CXX) doesn't print a "command not found" warning on every
# `make` invocation that doesn't even target lto-build.
lto-build:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	@CXX_VERSION="$$($(CXX) --version 2>&1 | head -1)"; \
	  case "$$CXX_VERSION" in *clang*) LTO_FLAG=-flto=thin ;; *) LTO_FLAG=-flto ;; esac; \
	  echo "LTO_FLAG=$$LTO_FLAG (cxx: $$CXX_VERSION, arch: $(LTO_ARCH))"; \
	  $(MAKE) arch=$(LTO_ARCH) EXE=bwa-mem3.lto EXTRA_CXXFLAGS="$(EXTRA_CXXFLAGS) $$LTO_FLAG -fno-semantic-interposition" CXX="$(CXX)" all
	@echo "LTO binary: bwa-mem3.lto (arch=$(LTO_ARCH))"
	# Drop variant-flagged objects from the shared cache so a subsequent
	# `make all` doesn't relink stale -flto / -fno-semantic-interposition
	# objects.
	rm -f src/*.o $(BWA_LIB)

lto-clean:
	rm -f bwa-mem3.lto

# Print the effective mimalloc setting. Used by CI and humans.
print-mimalloc-config:
	@echo "USE_MIMALLOC=$(USE_MIMALLOC)"

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -- src/*.cpp)

# DO NOT DELETE

src/FMI_search.o: src/bwa_madvise.h src/bwa_shm.h src/FMI_search.h
src/FMI_search.o: src/simd_compat.h src/read_index_ele.h src/utils.h
src/FMI_search.o: src/bntseq.h src/macro.h src/bwa.h src/bwt.h
src/FMI_search.o: src/memcpy_bwamem.h src/safestringlib.h
src/FMI_search.o: ext/safestringlib/include/safe_mem_lib.h
src/FMI_search.o: ext/safestringlib/include/safe_lib.h
src/FMI_search.o: ext/safestringlib/include/safe_types.h
src/FMI_search.o: ext/safestringlib/include/safe_lib_errno.h
src/FMI_search.o: ext/safestringlib/include/safe_str_lib.h src/profiling.h
src/FMI_search.o: src/libsais_build.h
src/bam_writer.o: ext/htslib/htslib/sam.h ext/htslib/htslib/hts.h
src/bam_writer.o: ext/htslib/htslib/hts_defs.h ext/htslib/htslib/hts_log.h
src/bam_writer.o: ext/htslib/htslib/kstring.h ext/htslib/htslib/kroundup.h
src/bam_writer.o: ext/htslib/htslib/hts_endian.h ext/htslib/htslib/kstring.h
src/bam_writer.o: src/bam_writer.h src/bwa.h src/bntseq.h src/bwt.h
src/bam_writer.o: src/macro.h src/bwamem.h src/kthread.h src/bandedSWA.h
src/bam_writer.o: src/simd_compat.h src/kernel_dispatch.h src/kstring.h
src/bam_writer.o: src/ksw.h src/kvec.h src/ksort.h src/utils.h
src/bam_writer.o: src/profiling.h src/FMI_search.h src/read_index_ele.h
src/bam_writer.o: src/cigar_util.h
src/bandedSWA.o: src/kernel_dispatch.h src/bandedSWA.h src/macro.h
src/bandedSWA.o: src/simd_compat.h
src/bntseq.o: src/bntseq.h src/utils.h src/macro.h src/kseq.h
src/bntseq.o: src/memcpy_bwamem.h src/safestringlib.h
src/bntseq.o: ext/safestringlib/include/safe_mem_lib.h
src/bntseq.o: ext/safestringlib/include/safe_lib.h
src/bntseq.o: ext/safestringlib/include/safe_types.h
src/bntseq.o: ext/safestringlib/include/safe_lib_errno.h
src/bntseq.o: ext/safestringlib/include/safe_str_lib.h src/khash.h
src/bwa.o: src/bntseq.h src/bwa.h src/bwt.h src/macro.h src/ksw.h
src/bwa.o: src/kernel_dispatch.h src/simd_compat.h src/utils.h src/kstring.h
src/bwa.o: src/kvec.h src/u8vec_scratch.h src/safestringlib.h
src/bwa.o: ext/safestringlib/include/safe_mem_lib.h
src/bwa.o: ext/safestringlib/include/safe_lib.h
src/bwa.o: ext/safestringlib/include/safe_types.h
src/bwa.o: ext/safestringlib/include/safe_lib_errno.h
src/bwa.o: ext/safestringlib/include/safe_str_lib.h src/kseq.h
src/bwa.o: src/memcpy_bwamem.h
src/bwa_shm.o: src/bwa_shm.h src/bwa.h src/bntseq.h src/bwt.h src/macro.h
src/bwa_shm.o: src/FMI_search.h src/simd_compat.h src/read_index_ele.h
src/bwa_shm.o: src/utils.h src/safestringlib.h
src/bwa_shm.o: ext/safestringlib/include/safe_mem_lib.h
src/bwa_shm.o: ext/safestringlib/include/safe_lib.h
src/bwa_shm.o: ext/safestringlib/include/safe_types.h
src/bwa_shm.o: ext/safestringlib/include/safe_lib_errno.h
src/bwa_shm.o: ext/safestringlib/include/safe_str_lib.h
src/bwamem.o: src/bwamem.h src/bwt.h src/bntseq.h src/bwa.h src/macro.h
src/bwamem.o: src/kthread.h src/bandedSWA.h src/simd_compat.h
src/bwamem.o: src/kernel_dispatch.h src/kstring.h src/ksw.h src/kvec.h
src/bwamem.o: src/ksort.h src/utils.h src/profiling.h src/FMI_search.h
src/bwamem.o: src/read_index_ele.h src/memcpy_bwamem.h src/safestringlib.h
src/bwamem.o: ext/safestringlib/include/safe_mem_lib.h
src/bwamem.o: ext/safestringlib/include/safe_lib.h
src/bwamem.o: ext/safestringlib/include/safe_types.h
src/bwamem.o: ext/safestringlib/include/safe_lib_errno.h
src/bwamem.o: ext/safestringlib/include/safe_str_lib.h src/bam_writer.h
src/bwamem.o: src/meth_bam.h src/u8vec_scratch.h src/sam_encode.h
src/bwamem.o: src/kbtree.h
src/bwamem_extra.o: src/bwa.h src/bntseq.h src/bwt.h src/macro.h src/bwamem.h
src/bwamem_extra.o: src/kthread.h src/bandedSWA.h src/simd_compat.h
src/bwamem_extra.o: src/kernel_dispatch.h src/kstring.h src/ksw.h src/kvec.h
src/bwamem_extra.o: src/ksort.h src/utils.h src/profiling.h src/FMI_search.h
src/bwamem_extra.o: src/read_index_ele.h
src/bwamem_pair.o: src/kstring.h src/bwamem.h src/bwt.h src/bntseq.h
src/bwamem_pair.o: src/bwa.h src/macro.h src/kthread.h src/bandedSWA.h
src/bwamem_pair.o: src/simd_compat.h src/kernel_dispatch.h src/ksw.h
src/bwamem_pair.o: src/kvec.h src/ksort.h src/utils.h src/profiling.h
src/bwamem_pair.o: src/FMI_search.h src/read_index_ele.h src/bam_writer.h
src/bwamem_pair.o: src/u8vec_scratch.h src/kswv.h
src/bwtindex.o: src/bntseq.h src/bwa.h src/bwt.h src/macro.h src/utils.h
src/bwtindex.o: src/FMI_search.h src/simd_compat.h src/read_index_ele.h
src/bwtindex.o: src/kseq.h src/memcpy_bwamem.h src/safestringlib.h
src/bwtindex.o: ext/safestringlib/include/safe_mem_lib.h
src/bwtindex.o: ext/safestringlib/include/safe_lib.h
src/bwtindex.o: ext/safestringlib/include/safe_types.h
src/bwtindex.o: ext/safestringlib/include/safe_lib_errno.h
src/bwtindex.o: ext/safestringlib/include/safe_str_lib.h src/system.h
src/fastmap.o: src/bwa_madvise.h src/fastmap.h src/bwa.h src/bntseq.h
src/fastmap.o: src/bwt.h src/macro.h src/bwamem.h src/kthread.h
src/fastmap.o: src/bandedSWA.h src/simd_compat.h src/kernel_dispatch.h
src/fastmap.o: src/kstring.h src/ksw.h src/kvec.h src/ksort.h src/utils.h
src/fastmap.o: src/profiling.h src/FMI_search.h src/read_index_ele.h
src/fastmap.o: src/kseq.h src/memcpy_bwamem.h src/safestringlib.h
src/fastmap.o: ext/safestringlib/include/safe_mem_lib.h
src/fastmap.o: ext/safestringlib/include/safe_lib.h
src/fastmap.o: ext/safestringlib/include/safe_types.h
src/fastmap.o: ext/safestringlib/include/safe_lib_errno.h
src/fastmap.o: ext/safestringlib/include/safe_str_lib.h src/bam_writer.h
src/fastmap.o: src/meth_bam.h src/meth_orig_ref.h src/bwa_shm.h
src/fm_index_writer.o: src/fm_index_writer.h src/FMI_search.h
src/fm_index_writer.o: src/simd_compat.h src/read_index_ele.h src/utils.h
src/fm_index_writer.o: src/bntseq.h src/macro.h src/bwa.h src/bwt.h
src/fm_index_writer.o: src/io_utils.h
src/index_prelude.o: src/index_prelude.h src/io_utils.h src/utils.h
src/index_prelude.o: src/packed_text.h
src/kopen.o: src/memcpy_bwamem.h src/safestringlib.h
src/kopen.o: ext/safestringlib/include/safe_mem_lib.h
src/kopen.o: ext/safestringlib/include/safe_lib.h
src/kopen.o: ext/safestringlib/include/safe_types.h
src/kopen.o: ext/safestringlib/include/safe_lib_errno.h
src/kopen.o: ext/safestringlib/include/safe_str_lib.h
src/kstring.o: src/kstring.h
src/ksw.o: src/kernel_dispatch.h src/simd_compat.h src/ksw.h src/macro.h
src/kswv.o: src/kernel_dispatch.h src/kswv.h src/macro.h src/ksw.h
src/kswv.o: src/simd_compat.h src/bandedSWA.h src/neon_utils.h
src/kthread.o: src/kthread.h src/macro.h src/bwamem.h src/bwt.h src/bntseq.h
src/kthread.o: src/bwa.h src/bandedSWA.h src/simd_compat.h
src/kthread.o: src/kernel_dispatch.h src/kstring.h src/ksw.h src/kvec.h
src/kthread.o: src/ksort.h src/utils.h src/profiling.h src/FMI_search.h
src/kthread.o: src/read_index_ele.h
src/libsais_build.o: src/libsais_build.h src/fm_index_writer.h
src/libsais_build.o: src/index_prelude.h src/io_utils.h src/utils.h
src/libsais_build.o: src/macro.h src/packed_text.h
src/libsais_build.o: ext/libsais/include/libsais.h
src/libsais_build.o: ext/libsais/include/libsais64.h
src/main.o: src/main.h src/kstring.h src/utils.h src/macro.h src/bandedSWA.h
src/main.o: src/simd_compat.h src/kernel_dispatch.h src/profiling.h
src/main.o: src/fastmap.h src/bwa.h src/bntseq.h src/bwt.h src/bwamem.h
src/main.o: src/kthread.h src/ksw.h src/kvec.h src/ksort.h src/FMI_search.h
src/main.o: src/read_index_ele.h src/kseq.h src/memcpy_bwamem.h
src/main.o: src/safestringlib.h ext/safestringlib/include/safe_mem_lib.h
src/main.o: ext/safestringlib/include/safe_lib.h
src/main.o: ext/safestringlib/include/safe_types.h
src/main.o: ext/safestringlib/include/safe_lib_errno.h
src/main.o: ext/safestringlib/include/safe_str_lib.h src/simd_dispatch.h
src/main.o: src/version.h src/bwa_shm.h ext/mimalloc/include/mimalloc.h
src/memcpy_bwamem.o: src/memcpy_bwamem.h src/safestringlib.h
src/memcpy_bwamem.o: ext/safestringlib/include/safe_mem_lib.h
src/memcpy_bwamem.o: ext/safestringlib/include/safe_lib.h
src/memcpy_bwamem.o: ext/safestringlib/include/safe_types.h
src/memcpy_bwamem.o: ext/safestringlib/include/safe_lib_errno.h
src/memcpy_bwamem.o: ext/safestringlib/include/safe_str_lib.h
src/meth_bam.o: ext/htslib/htslib/sam.h ext/htslib/htslib/hts.h
src/meth_bam.o: ext/htslib/htslib/hts_defs.h ext/htslib/htslib/hts_log.h
src/meth_bam.o: ext/htslib/htslib/kstring.h ext/htslib/htslib/kroundup.h
src/meth_bam.o: ext/htslib/htslib/hts_endian.h ext/htslib/htslib/kstring.h
src/meth_bam.o: src/meth_bam.h src/bwa.h src/bntseq.h src/bwt.h src/macro.h
src/meth_bam.o: src/bwamem.h src/kthread.h src/bandedSWA.h src/simd_compat.h
src/meth_bam.o: src/kernel_dispatch.h src/kstring.h src/ksw.h src/kvec.h
src/meth_bam.o: src/ksort.h src/utils.h src/profiling.h src/FMI_search.h
src/meth_bam.o: src/read_index_ele.h src/bam_writer.h src/cigar_util.h
src/meth_bam.o: src/meth_orig_ref.h src/meth_xm.h src/version.h
src/meth_orig_ref.o: src/meth_orig_ref.h src/bntseq.h src/meth_bam.h
src/meth_orig_ref.o: src/bwa.h src/bwt.h src/macro.h src/bwamem.h
src/meth_orig_ref.o: src/kthread.h src/bandedSWA.h src/simd_compat.h
src/meth_orig_ref.o: src/kernel_dispatch.h src/kstring.h src/ksw.h src/kvec.h
src/meth_orig_ref.o: src/ksort.h src/utils.h src/profiling.h src/FMI_search.h
src/meth_orig_ref.o: src/read_index_ele.h
src/meth_xm.o: src/meth_xm.h src/meth_orig_ref.h src/bntseq.h src/meth_bam.h
src/meth_xm.o: src/bwa.h src/bwt.h src/macro.h src/bwamem.h src/kthread.h
src/meth_xm.o: src/bandedSWA.h src/simd_compat.h src/kernel_dispatch.h
src/meth_xm.o: src/kstring.h src/ksw.h src/kvec.h src/ksort.h src/utils.h
src/meth_xm.o: src/profiling.h src/FMI_search.h src/read_index_ele.h
src/packed_text.o: src/packed_text.h src/utils.h
src/profiling.o: src/macro.h src/bwa.h src/bntseq.h src/bwt.h src/profiling.h
src/read_index_ele.o: src/read_index_ele.h src/utils.h src/bntseq.h
src/read_index_ele.o: src/macro.h src/safestringlib.h
src/read_index_ele.o: ext/safestringlib/include/safe_mem_lib.h
src/read_index_ele.o: ext/safestringlib/include/safe_lib.h
src/read_index_ele.o: ext/safestringlib/include/safe_types.h
src/read_index_ele.o: ext/safestringlib/include/safe_lib_errno.h
src/read_index_ele.o: ext/safestringlib/include/safe_str_lib.h
src/read_index_ele.o: src/bwa_madvise.h src/bwa_shm.h
src/sam_encode.o: src/sam_encode.h src/kernel_dispatch.h
src/simd_dispatch.o: src/simd_dispatch.h src/bandedSWA.h src/macro.h
src/simd_dispatch.o: src/simd_compat.h src/kernel_dispatch.h src/kswv.h
src/simd_dispatch.o: src/ksw.h src/sam_encode.h
src/system.o: src/system.h
src/utils.o: src/utils.h src/ksort.h src/kseq.h src/memcpy_bwamem.h
src/utils.o: src/safestringlib.h ext/safestringlib/include/safe_mem_lib.h
src/utils.o: ext/safestringlib/include/safe_lib.h
src/utils.o: ext/safestringlib/include/safe_types.h
src/utils.o: ext/safestringlib/include/safe_lib_errno.h
src/utils.o: ext/safestringlib/include/safe_str_lib.h
