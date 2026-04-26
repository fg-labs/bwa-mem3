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

EXE=		bwa-mem2
#CXX=		icpc
ifeq ($(CXX), icpc)
	CC= icc
else ifeq ($(CXX), g++)
	CC=gcc
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
# Override with USE_MIMALLOC=0 to build a stock bwa-mem2 without mimalloc.
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

MEM_FLAGS=	-DSAIS=1
CPPFLAGS+=	-DENABLE_PREFETCH -DV17=1 -DMATE_SORT=0 $(MEM_FLAGS)

# Version string for `bwa-mem2 version` and the @PG VN: field. Prefer
# `git describe` (e.g. v2.3-30-g61813ef, with -dirty suffix for modified
# trees) so the stamped version always reflects the actual build. Fall
# back to a static tag for source-tarball / shallow-clone builds.
FG_LABS_VERSION_FALLBACK := 2.3-fg-labs
VERSION_STRING := $(shell git describe --tags --dirty 2>/dev/null || echo $(FG_LABS_VERSION_FALLBACK))
INCLUDES+=   -Isrc -Iext/safestringlib/include -Iext/htslib
ifeq ($(USE_MIMALLOC),1)
    INCLUDES += -Iext/mimalloc/include
endif
LIBS=		-lpthread -lm -lz -L. -lbwa -Lext/safestringlib -lsafestring -Lext/htslib -lhts $(STATIC_GCC) $(LIBS_EXTRA)
OBJS=		src/fastmap.o src/bwtindex.o src/utils.o src/memcpy_bwamem.o src/kthread.o \
			src/kstring.o src/ksw.o src/bntseq.o src/bwamem.o src/profiling.o src/bandedSWA.o \
			src/FMI_search.o src/read_index_ele.o src/bwamem_pair.o src/kswv.o src/bwa.o \
			src/bwamem_extra.o src/kopen.o src/bam_writer.o src/meth_bam.o
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

CXXFLAGS+=	-g -O3 -std=gnu++14 -fpermissive $(ARCH_FLAGS) $(ASAN_FLAGS) $(EXTRA_CXXFLAGS) #-Wall ##-xSSE2

# Control build flag for the batched mate-rescue SW port on ARM.
# When set (e.g. `make arm64 DISABLE_BATCHED_MATESW=1`), the source gate for
# the new batched path falls through to the legacy scalar mem_sam_pe. Used by
# the proto-neon-kswv CI to A/B the same commit with the port on vs. off.
# Pass the caller-supplied value through verbatim so `DISABLE_BATCHED_MATESW=0`
# still selects the batched path (ifdef would be true even for =0).
ifneq ($(strip $(DISABLE_BATCHED_MATESW)),)
    CPPFLAGS += -DDISABLE_BATCHED_MATESW=$(DISABLE_BATCHED_MATESW)
endif

.PHONY:all clean depend multi print-mimalloc-config kswv_selftest kswv_nrow_zero_test test FORCE
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
	@printf '#ifndef BWA_MEM2_VERSION_H\n#define BWA_MEM2_VERSION_H\n#define PACKAGE_VERSION "%s"\n#endif\n' '$(VERSION_STRING)' > $@.tmp
	@if ! cmp -s $@.tmp $@ 2>/dev/null; then mv $@.tmp $@; else rm -f $@.tmp; fi

src/main.o: src/version.h

multi: $(if $(filter 1,$(USE_MIMALLOC)),$(MIMALLOC_LIB))
ifneq ($(IS_ARM),)
	@echo "ARM64 detected - building single arm64 binary instead of multi"
	$(MAKE) arm64
else
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse41    EXE=bwa-mem2.sse41    CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse42    EXE=bwa-mem2.sse42    CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx    EXE=bwa-mem2.avx    CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx2   EXE=bwa-mem2.avx2     CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx512bw EXE=bwa-mem2.avx512bw CXX=$(CXX) all
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -Wall -O3 src/runsimd.cpp -Iext/safestringlib/include -Lext/safestringlib/ -lsafestring $(STATIC_GCC) -o bwa-mem2
endif

# ARM64/Apple Silicon build target - single binary, no multi-binary launcher needed
arm64:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=arm64 EXE=bwa-mem2.arm64 CXX=$(CXX) all
	ln -sf bwa-mem2.arm64 bwa-mem2


$(EXE):$(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) $(if $(filter 1,$(USE_MIMALLOC)),$(MIMALLOC_LIB)) src/main.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) src/main.o $(BWA_LIB) $(LIBS) $(MIMALLOC_LDFLAGS) -o $@

# kswv self-consistency test: batched SIMD kswv vs scalar ksw_align2 reference.
# Built by the proto-neon-kswv CI workflow; runnable standalone.
kswv_selftest: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) test/kswv_selftest.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) test/kswv_selftest.o $(BWA_LIB) $(LIBS) -o $@

# Regression test for issue 38 / upstream PR 289: exercises an all-len1==0
# batch that drives each SIMD kswv kernel through the nrow==0 path. Without
# the post-loop `if (i > 0)` guard, the rowMax store writes SIMD_WIDTH* bytes
# before the allocation and aborts at a later allocator operation; under
# asan the write is reported directly.
kswv_nrow_zero_test: $(BWA_LIB) $(SAFE_STR_LIB) $(HTS_LIB) test/kswv_nrow_zero_test.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) test/kswv_nrow_zero_test.o $(BWA_LIB) $(LIBS) -o $@

# Run the in-tree tests. Currently kswv_selftest + kswv_nrow_zero_test;
# extend as more land.
test: kswv_selftest kswv_nrow_zero_test
	./kswv_selftest
	./kswv_nrow_zero_test

test/kswv_selftest.o: test/kswv_selftest.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

test/kswv_nrow_zero_test.o: test/kswv_nrow_zero_test.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

$(BWA_LIB):$(OBJS)
	ar rcs $(BWA_LIB) $(OBJS)

$(HTS_LIB):
	cd ext/htslib && \
	    ([ -f Makefile ] || (autoreconf -i && \
	        ./configure --disable-lzma --disable-libcurl --disable-gcs \
	                    --disable-s3 --disable-plugins --disable-bz2)) && \
	    $(MAKE) libhts.a

# On macOS, safestringlib needs stdlib.h for abort() and the memset_s
# declaration conflicts with macOS C11 Annex K (different signature).
SAFE_EXTRA_CFLAGS =
ifeq ($(UNAME_S),Darwin)
    SAFE_EXTRA_CFLAGS = -include stdlib.h -Dmemset_s=_safestringlib_memset_s
endif

$(SAFE_STR_LIB):
	cd ext/safestringlib/ && $(MAKE) clean && $(MAKE) CC=$(CC) CFLAGS="-Iinclude -Isafeclib $(SAFE_EXTRA_CFLAGS) -fstack-protector-strong -fPIE -fPIC -O2" directories libsafestring.a

# htslib: minimal configure (no lzma/bz2/curl/S3/GCS/plugins), zlib only.
# Guard on config.mk (only created by ./configure) rather than Makefile, which
# is checked into the htslib tree and would make the guard a no-op.
$(HTS_LIB):
	cd ext/htslib && \
	    ([ -f config.mk ] || (autoreconf -i && \
	        ./configure --disable-lzma --disable-libcurl --disable-gcs \
	                    --disable-s3 --disable-plugins --disable-bz2)) && \
	    $(MAKE) libhts.a

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

clean:
	rm -fr src/*.o src/version.h test/*.o $(BWA_LIB) $(EXE) kswv_selftest kswv_nrow_zero_test bwa-mem2.sse41 bwa-mem2.sse42 bwa-mem2.avx bwa-mem2.avx2 bwa-mem2.avx512bw bwa-mem2.arm64 bwa-mem2.pgo bwa-mem2.pgo-instr bwa-mem2.pgo.* bwa-mem2.pgo-instr.*
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
#   PGO_ARCH=arm64  -> bwa-mem2.pgo-instr,    bwa-mem2.pgo
#   PGO_ARCH=avx2   -> bwa-mem2.pgo-instr.avx2, bwa-mem2.pgo.avx2
ifneq ($(IS_ARM),)
    PGO_ARCH ?= arm64
else
    PGO_ARCH ?= native
endif
PGO_PROFILE_DIR ?= pgo_profiles

# Output names: keep the bare names when PGO_ARCH is the default arm64
# (backward-compat); arch-suffix otherwise so per-arch outputs don't collide.
ifeq ($(PGO_ARCH),arm64)
    PGO_INSTR_EXE = bwa-mem2.pgo-instr
    PGO_FINAL_EXE = bwa-mem2.pgo
else
    PGO_INSTR_EXE = bwa-mem2.pgo-instr.$(PGO_ARCH)
    PGO_FINAL_EXE = bwa-mem2.pgo.$(PGO_ARCH)
endif

pgo-generate:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=$(PGO_ARCH) EXE=$(PGO_INSTR_EXE) EXTRA_CXXFLAGS="-fprofile-generate=$(PGO_PROFILE_DIR)" CXX=$(CXX) all
	@echo "PGO instrumented binary built: $(PGO_INSTR_EXE) (arch=$(PGO_ARCH), profile dir=$(PGO_PROFILE_DIR))"
	@echo "Run training workload with $(PGO_INSTR_EXE), then: make pgo-use PGO_ARCH=$(PGO_ARCH) PGO_PROFILE_DIR=$(PGO_PROFILE_DIR)"

pgo-use:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=$(PGO_ARCH) EXE=$(PGO_FINAL_EXE) EXTRA_CXXFLAGS="-fprofile-use=$(PGO_PROFILE_DIR) -fprofile-correction" CXX=$(CXX) all
	@echo "PGO optimized binary built: $(PGO_FINAL_EXE) (arch=$(PGO_ARCH))"

pgo-clean:
	rm -rf $(PGO_PROFILE_DIR) bwa-mem2.pgo-instr bwa-mem2.pgo bwa-mem2.pgo-instr.* bwa-mem2.pgo.*

# Print the effective mimalloc setting. Used by CI and humans.
print-mimalloc-config:
	@echo "USE_MIMALLOC=$(USE_MIMALLOC)"

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(CPPFLAGS) -I. -- src/*.cpp)

# DO NOT DELETE

src/FMI_search.o: src/FMI_search.h src/bntseq.h src/read_index_ele.h
src/FMI_search.o: src/utils.h src/macro.h src/bwa.h src/bwt.h src/sais.h
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
src/meth_bam.o: src/meth_bam.h src/bwamem.h src/bwa.h src/bntseq.h
