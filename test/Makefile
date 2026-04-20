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


EXE=		fmi_test smem2_test bwt_seed_strategy_test sa2ref_test xeonbsw

# Match the parent Makefile's compiler + flag selection. Detect ARM
# (Apple Silicon arm64 and Linux aarch64) and wire in OpenMP, sse2neon
# SSE-to-NEON translation, and (on macOS) Accelerate. On Linux x86 the
# original Intel toolchain used icpc; this now defaults to g++ so CI
# (which has no icpc) can build. Override with e.g. `make CXX=icpc` or
# `make CXX=clang++` to restore or change.
UNAME_M := $(shell uname -m)
UNAME_S := $(shell uname -s)
IS_ARM  := $(filter $(UNAME_M),arm64 aarch64)

ifneq ($(IS_ARM),)
    # ARM build: requires libomp for OpenMP (brew install libomp on macOS;
    # libgomp via system compiler on Linux aarch64) and sse2neon to
    # translate the SSE intrinsics used through ../src/simd_compat.h.
    CXX ?= clang++
    ifeq ($(UNAME_S),Darwin)
        LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
        ifeq ($(strip $(LIBOMP_PREFIX)),)
            $(error libomp not found: install with `brew install libomp`, or set LIBOMP_PREFIX to its install prefix)
        endif
        CXXFLAGS ?= -std=c++11 -O2 -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
        LIBS ?= -L.. -L$(LIBOMP_PREFIX)/lib -lomp -lbwa -lz -framework Accelerate ../ext/safestringlib/libsafestring.a
    else
        CXXFLAGS ?= -std=c++11 -O2 -fopenmp
        LIBS ?= -L.. -fopenmp -lbwa -lz ../ext/safestringlib/libsafestring.a
    endif
else
    CXX ?=      g++
    CXXFLAGS ?= -std=c++11 -fopenmp -mtune=native -march=native
    LIBS ?=     -L.. -fopenmp -lbwa -lz ../ext/safestringlib/libsafestring.a
endif

# CPPFLAGS / INCLUDES: add sse2neon SSE-emulation macros + include path
# on ARM so ../src/simd_compat.h can translate the SSE intrinsics used
# by fmi_test.cpp et al. Add APPLE_SILICON + CACHE_LINE_BYTES=128 to
# match the parent Makefile's ARM flags.
ifneq ($(IS_ARM),)
    CPPFLAGS = -DENABLE_PREFETCH -D__SSE__=1 -D__SSE2__=1 -D__SSE3__=1 -D__SSSE3__=1 -D__SSE4_1__=1 -D__SSE4_2__=1 -DAPPLE_SILICON=1 -DCACHE_LINE_BYTES=128
    INCLUDES = -I../src -I../ext/safestringlib/include -I../ext/sse2neon
else
    CPPFLAGS = -DENABLE_PREFETCH
    INCLUDES = -I../src -I../ext/safestringlib/include
endif

.PHONY:all clean depend
.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

fmi_test:fmi_test.o
	$(CXX) -o $@ $^ $(LIBS)

smem2_test: smem2_test.o
	$(CXX) -o $@ $^ $(LIBS)

sa2ref_test: sa2ref_test.o
	$(CXX) -o $@ $^ $(LIBS)

bwt_seed_strategy_test:bwt_seed_strategy_test.o
	$(CXX) -o $@ $^ $(LIBS)

xeonbsw:main_banded.o
	$(CXX) -o $@ $^ $(LIBS)

clean:
	rm -fr *.o $(EXE)


# DO NOT DELETE

bwt_seed_strategy_test.o: ../src/FMI_search.h ../src/bntseq.h ../src/read_index_ele.h
bwt_seed_strategy_test.o: ../src/bwa.h ../src/bwt.h ../src/utils.h ../src/macro.h
fmi_test.o: ../src/FMI_search.h ../src/bntseq.h ../src/read_index_ele.h
fmi_test.o: ../src/bwa.h ../src/bwt.h ../src/utils.h ../src/macro.h
main_banded.o: ../src/bandedSWA.h ../src/macro.h
sa2ref_test.o: ../src/FMI_search.h ../src/bntseq.h ../src/read_index_ele.h
sa2ref_test.o: ../src/bwa.h ../src/bwt.h ../src/utils.h ../src/macro.h
smem2_test.o: ../src/FMI_search.h ../src/bntseq.h ../src/read_index_ele.h
smem2_test.o: ../src/bwa.h ../src/bwt.h ../src/utils.h ../src/macro.h
