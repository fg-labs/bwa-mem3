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
FG_LABS_VERSION_FALLBACK := 0.1.0-pre
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
OBJS=		src/fastmap.o src/bwtindex.o src/utils.o src/memcpy_bwamem.o src/kthread.o \
			src/kstring.o src/ksw.o src/bntseq.o src/bwamem.o src/profiling.o src/bandedSWA.o \
			src/FMI_search.o src/read_index_ele.o src/bwamem_pair.o src/kswv.o src/bwa.o \
			src/bwamem_extra.o src/kopen.o src/bam_writer.o src/meth_bam.o \
			src/packed_text.o src/fm_index_writer.o src/index_prelude.o \
			src/system.o src/libsais_build.o \
			src/bwa_shm.o
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
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-xCORE-AVX512
	else
		ARCH_FLAGS=-mavx512f -mavx512bw
	endif
else ifeq ($(arch),avx512bw)
	# Explicit BW target: double the lane width vs AVX2 (64x8-bit / 32x16-bit).
	# AVX-512BW implies AVX-512F; -mavx512bw alone enables BW+F on gcc/clang
	# but we list both flags for clarity.
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-xCORE-AVX512
	else
		ARCH_FLAGS=-mavx512f -mavx512bw
	endif
else ifeq ($(arch),native)
	ARCH_FLAGS=-march=native
else ifneq ($(arch),)
# To provide a different architecture flag like -march=core-avx2.
	ARCH_FLAGS=$(arch)
else
myall:multi
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

.PHONY:all clean depend multi print-mimalloc-config kswv_nrow_zero_test shm_section_find_test shm_pack_round_trip_test test FORCE pgo-generate pgo-use pgo-clean profile-build profile-clean lto-build lto-clean
.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

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

multi: $(if $(filter 1,$(USE_MIMALLOC)),$(MIMALLOC_LIB))
ifneq ($(IS_ARM),)
	@echo "ARM64 detected - building single arm64 binary instead of multi"
	$(MAKE) arm64
else
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse41    EXE=bwa-mem3.sse41    CXX="$(CXX)" all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse42    EXE=bwa-mem3.sse42    CXX="$(CXX)" all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx    EXE=bwa-mem3.avx    CXX="$(CXX)" all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx2   EXE=bwa-mem3.avx2     CXX="$(CXX)" all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx512bw EXE=bwa-mem3.avx512bw CXX="$(CXX)" all
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -Wall -O3 src/runsimd.cpp -Iext/safestringlib/include -Lext/safestringlib/ -lsafestring $(STATIC_GCC) -o bwa-mem3
endif

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
kswv_nrow_zero_test: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) test/kswv_nrow_zero_test.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) test/kswv_nrow_zero_test.o $(BWA_LIB) $(LIBS) -o $@

# Build the test binaries with the same ARCH_FLAGS as libbwa.a so the
# test binary's kswv.h preprocessor state (SIMD_WIDTH8, BWA_TESTS_HAVE_KSWV)
# matches what libbwa.a was compiled with. Consumed by ci.yml so that e.g.
# the sse41 matrix row builds test/framework with -msse4.1 only (matching
# libbwa.a, which then lacks kswv::getScores8 — the BWA_TESTS_HAVE_KSWV
# macro guards the test away).
.PHONY: test-binaries
# $(SAFE_STR_LIB) is a real link-time dep: test/Makefile's
# bwa_mem3_tests_unit recipe references ../ext/safestringlib/libsafestring.a
# directly. Without this prereq, callers that skip the bwa-mem3 binary
# build (which builds it as a side-effect of $(EXE) deps) link-fail.
test-binaries: $(BWA_LIB) $(SAFE_STR_LIB)
	$(MAKE) -C test framework unit integration \
	    CXX="$(CXX)" \
	    COVERAGE=$(COVERAGE) \
	    ARCH_FLAGS_FROM_PARENT='$(ARCH_FLAGS)'

shm_section_find_test: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) $(LIBSAIS_OBJS) test/shm_section_find_test.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) test/shm_section_find_test.o $(BWA_LIB) $(LIBSAIS_OBJS) $(LIBS) -o $@

shm_pack_round_trip_test: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) $(LIBSAIS_OBJS) test/shm_pack_round_trip_test.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) test/shm_pack_round_trip_test.o $(BWA_LIB) $(LIBSAIS_OBJS) $(LIBS) -o $@

test/shm_pack_round_trip_test.o: test/shm_pack_round_trip_test.cpp

# Run the in-tree tests via the unit-test harness in test/, plus the
# standalone regressions (kswv_nrow_zero_test + shm_section_find_test).
# shm_pack_round_trip_test runs via test/shm_pack_round_trip_test.sh which
# builds the phiX index first; invoked from test/run_unit_tests.sh.
test: test-binaries kswv_nrow_zero_test shm_section_find_test
	./test/bwa_mem3_tests_unit
	./test/bwa_mem3_tests_integration
	./kswv_nrow_zero_test
	./shm_section_find_test

test/kswv_nrow_zero_test.o: test/kswv_nrow_zero_test.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

test/shm_section_find_test.o: test/shm_section_find_test.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

$(BWA_LIB):$(OBJS)
	ar rcs $(BWA_LIB) $(OBJS)

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
	rm -fr src/*.o src/version.h test/*.o $(BWA_LIB) $(EXE) kswv_nrow_zero_test shm_section_find_test shm_pack_round_trip_test bwa-mem3.sse41 bwa-mem3.sse42 bwa-mem3.avx bwa-mem3.avx2 bwa-mem3.avx512bw bwa-mem3.arm64
	rm -f $(LIBSAIS_OBJS)
	rm -f src/*.gcno src/*.gcda
	$(MAKE) -C test clean
	cd ext/safestringlib/ && $(MAKE) clean
	-[ -f ext/htslib/config.mk ] && cd ext/htslib && $(MAKE) distclean
	rm -rf $(MIMALLOC_BUILD)

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

src/FMI_search.o: src/FMI_search.h src/bntseq.h src/read_index_ele.h
src/FMI_search.o: src/utils.h src/macro.h src/bwa.h src/bwt.h
src/FMI_search.o: src/libsais_build.h
src/libsais_build.o: src/libsais_build.h src/fm_index_writer.h src/index_prelude.h
src/libsais_build.o: src/packed_text.h src/utils.h src/macro.h
src/libsais_build.o: ext/libsais/include/libsais.h ext/libsais/include/libsais64.h
src/bandedSWA.o: src/bandedSWA.h src/macro.h
src/bntseq.o: src/bntseq.h src/utils.h src/macro.h src/kseq.h src/khash.h
src/bwa.o: src/bntseq.h src/bwa.h src/bwt.h src/macro.h src/ksw.h src/utils.h
src/bwa.o: src/kstring.h src/kvec.h src/kseq.h
src/bwamem.o: src/bwamem.h src/bwt.h src/bntseq.h src/bwa.h src/macro.h
src/bwamem.o: src/kthread.h src/bandedSWA.h src/kstring.h src/ksw.h
src/bwamem.o: src/kvec.h src/ksort.h src/utils.h src/profiling.h
src/bwamem.o: src/FMI_search.h src/read_index_ele.h src/kbtree.h
src/bwamem_extra.o: src/bwa.h src/bntseq.h src/bwt.h src/macro.h src/bwamem.h
src/bwamem_extra.o: src/kthread.h src/bandedSWA.h src/kstring.h src/ksw.h
src/bwamem_extra.o: src/kvec.h src/ksort.h src/utils.h src/profiling.h
src/bwamem_extra.o: src/FMI_search.h src/read_index_ele.h
src/bwamem_pair.o: src/kstring.h src/bwamem.h src/bwt.h src/bntseq.h
src/bwamem_pair.o: src/bwa.h src/macro.h src/kthread.h src/bandedSWA.h
src/bwamem_pair.o: src/ksw.h src/kvec.h src/ksort.h src/utils.h
src/bwamem_pair.o: src/profiling.h src/FMI_search.h src/read_index_ele.h
src/bwamem_pair.o: src/kswv.h
src/bwtindex.o: src/bntseq.h src/bwa.h src/bwt.h src/macro.h src/utils.h
src/bwtindex.o: src/FMI_search.h src/read_index_ele.h
src/fastmap.o: src/fastmap.h src/bwa.h src/bntseq.h src/bwt.h src/macro.h
src/fastmap.o: src/bwamem.h src/kthread.h src/bandedSWA.h src/kstring.h
src/fastmap.o: src/ksw.h src/kvec.h src/ksort.h src/utils.h src/profiling.h
src/fastmap.o: src/FMI_search.h src/read_index_ele.h src/kseq.h
src/kstring.o: src/kstring.h
src/ksw.o: src/ksw.h src/macro.h
src/kswv.o: src/kswv.h src/macro.h src/ksw.h src/bandedSWA.h
src/kthread.o: src/kthread.h src/macro.h src/bwamem.h src/bwt.h src/bntseq.h
src/kthread.o: src/bwa.h src/bandedSWA.h src/kstring.h src/ksw.h src/kvec.h
src/kthread.o: src/ksort.h src/utils.h src/profiling.h src/FMI_search.h
src/kthread.o: src/read_index_ele.h
src/main.o: src/main.h src/kstring.h src/utils.h src/macro.h src/bandedSWA.h
src/main.o: src/profiling.h
src/profiling.o: src/macro.h
src/read_index_ele.o: src/read_index_ele.h src/utils.h src/bntseq.h
src/read_index_ele.o: src/macro.h
src/utils.o: src/utils.h src/ksort.h src/kseq.h
src/memcpy_bwamem.o: src/memcpy_bwamem.h
src/bam_writer.o: src/bam_writer.h src/bwamem.h src/bwa.h src/bntseq.h
src/meth_bam.o: src/meth_bam.h src/bwamem.h src/bwa.h src/bntseq.h src/version.h
