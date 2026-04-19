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

# Detect architecture
UNAME_M := $(shell uname -m)
UNAME_S := $(shell uname -s)
# Treat macOS ("arm64") and Linux ("aarch64") as the same ARM build target.
IS_ARM := $(filter $(UNAME_M),arm64 aarch64)

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
INCLUDES+=   -Isrc -Iext/safestringlib/include
LIBS=		-lpthread -lm -lz -L. -lbwa -Lext/safestringlib -lsafestring $(STATIC_GCC) $(LIBS_EXTRA)
OBJS=		src/fastmap.o src/bwtindex.o src/utils.o src/memcpy_bwamem.o src/kthread.o \
			src/kstring.o src/ksw.o src/bntseq.o src/bwamem.o src/profiling.o src/bandedSWA.o \
			src/FMI_search.o src/read_index_ele.o src/bwamem_pair.o src/kswv.o src/bwa.o \
			src/bwamem_extra.o src/kopen.o
BWA_LIB=    libbwa.a
SAFE_STR_LIB=    ext/safestringlib/libsafestring.a

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
		ARCH_FLAGS=-mavx512bw
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

CXXFLAGS+=	-g -O3 -fpermissive $(ARCH_FLAGS) #-Wall ##-xSSE2

.PHONY:all clean depend multi
.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

multi:
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
	$(MAKE) arch=avx512 EXE=bwa-mem2.avx512bw CXX=$(CXX) all
	$(CXX) -Wall -O3 src/runsimd.cpp -Iext/safestringlib/include -Lext/safestringlib/ -lsafestring $(STATIC_GCC) -o bwa-mem2
endif

# ARM64/Apple Silicon build target - single binary, no multi-binary launcher needed
arm64:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=arm64 EXE=bwa-mem2.arm64 CXX=$(CXX) all
	ln -sf bwa-mem2.arm64 bwa-mem2


$(EXE):$(BWA_LIB) $(SAFE_STR_LIB) src/main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/main.o $(BWA_LIB) $(LIBS) -o $@

$(BWA_LIB):$(OBJS)
	ar rcs $(BWA_LIB) $(OBJS)

# On macOS, safestringlib needs stdlib.h for abort() and the memset_s
# declaration conflicts with macOS C11 Annex K (different signature).
SAFE_EXTRA_CFLAGS =
ifeq ($(UNAME_S),Darwin)
    SAFE_EXTRA_CFLAGS = -include stdlib.h -Dmemset_s=_safestringlib_memset_s
endif

$(SAFE_STR_LIB):
	cd ext/safestringlib/ && $(MAKE) clean && $(MAKE) CC=$(CC) CFLAGS="-Iinclude -Isafeclib $(SAFE_EXTRA_CFLAGS) -fstack-protector-strong -fPIE -fPIC -O2" directories libsafestring.a

clean:
	rm -fr src/*.o $(BWA_LIB) $(EXE) bwa-mem2.sse41 bwa-mem2.sse42 bwa-mem2.avx bwa-mem2.avx2 bwa-mem2.avx512bw bwa-mem2.arm64
	cd ext/safestringlib/ && $(MAKE) clean

# Profile-Guided Optimization (PGO) targets for Apple Silicon
# Usage: make pgo-generate && <run training workload> && make pgo-use
PGO_PROFILE_DIR=pgo_profiles

pgo-generate:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=arm64 EXE=bwa-mem2.pgo-instr CXXFLAGS="-g -O3 -fpermissive $(ARCH_FLAGS) -fprofile-generate=$(PGO_PROFILE_DIR)" CXX=$(CXX) all
	@echo "PGO instrumented binary built. Run training workload with bwa-mem2.pgo-instr"
	@echo "Then run: make pgo-use"

pgo-use:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=arm64 EXE=bwa-mem2.pgo CXXFLAGS="-g -O3 -fpermissive $(ARCH_FLAGS) -fprofile-use=$(PGO_PROFILE_DIR) -fprofile-correction" CXX=$(CXX) all
	@echo "PGO optimized binary built: bwa-mem2.pgo"

pgo-clean:
	rm -rf $(PGO_PROFILE_DIR) bwa-mem2.pgo-instr bwa-mem2.pgo

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
