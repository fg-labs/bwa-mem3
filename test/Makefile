##*****************************************************************************
##  bwa-mem2 test build
##
##  This Makefile produces:
##    - test/framework/libbwa_test_framework.a    (shared test-helper lib)
##    - test/bwa_mem2_tests_unit                  (doctest unit binary)
##    - test/bwa_mem2_tests_integration           (doctest integration binary)
##    - test/fmi_test, smem2_test, ...            (legacy per-binary tests,
##                                                  migrated one-by-one in
##                                                  follow-up PRs)
##
##  Arch flags (CXX, CXXFLAGS, CPPFLAGS, INCLUDES, LIBS) inherit from the
##  parent Makefile when invoked via `make -C test`. Override on the command
##  line for local builds that don't go through the parent (e.g. macOS ARM
##  direct builds — see the ARM branch below).
##*****************************************************************************

EXE = fmi_test smem2_test bwt_seed_strategy_test sa2ref_test xeonbsw smem_lockstep_parity_test header_insert_test packed_text_test system_test

# Doctest binaries — rules populated in Task 14 (unit) and Task 15
# (integration). Listed separately from $(EXE) so `all` and
# `run_unit_tests.sh` keep working during the framework-scaffolding PR.
DOCTEST_BINS = bwa_mem2_tests_unit bwa_mem2_tests_integration

UNAME_M := $(shell uname -m)
UNAME_S := $(shell uname -s)
IS_ARM  := $(filter $(UNAME_M),arm64 aarch64)

ifneq ($(IS_ARM),)
    CXX ?= clang++
    ifeq ($(UNAME_S),Darwin)
        LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
        ifeq ($(strip $(LIBOMP_PREFIX)),)
            $(error libomp not found: install with `brew install libomp`, or set LIBOMP_PREFIX to its install prefix)
        endif
        CXXFLAGS ?= -std=c++11 -O2 -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
        LIBS     ?= -L.. -L$(LIBOMP_PREFIX)/lib -lomp -lbwa ../ext/libsais/src/libsais.o ../ext/libsais/src/libsais64.o -lz -framework Accelerate ../ext/safestringlib/libsafestring.a
    else
        CXXFLAGS ?= -std=c++11 -O2 -fopenmp
        LIBS     ?= -L.. -fopenmp -lbwa ../ext/libsais/src/libsais.o ../ext/libsais/src/libsais64.o -lz ../ext/safestringlib/libsafestring.a
    endif
else
    CXX      ?= g++
    CXXFLAGS ?= -std=c++11 -fopenmp -mtune=native -march=native
    LIBS     ?= -L.. -fopenmp -lbwa ../ext/libsais/src/libsais.o ../ext/libsais/src/libsais64.o -lz ../ext/safestringlib/libsafestring.a
endif

ifneq ($(IS_ARM),)
    CPPFLAGS = -DENABLE_PREFETCH -D__SSE__=1 -D__SSE2__=1 -D__SSE3__=1 -D__SSSE3__=1 -D__SSE4_1__=1 -D__SSE4_2__=1 -DAPPLE_SILICON=1 -DCACHE_LINE_BYTES=128
    INCLUDES = -I../src -I../ext/safestringlib/include -I../ext/sse2neon -I../ext
else
    CPPFLAGS = -DENABLE_PREFETCH
    INCLUDES = -I../src -I../ext/safestringlib/include -I../ext
endif

FRAMEWORK_DIR  := framework
FRAMEWORK_SRCS := $(FRAMEWORK_DIR)/scoring.cpp \
                  $(FRAMEWORK_DIR)/seqpair_gen.cpp \
                  $(FRAMEWORK_DIR)/seqpair_batch.cpp \
                  $(FRAMEWORK_DIR)/ksw_runner.cpp \
                  $(FRAMEWORK_DIR)/kswv_runner.cpp \
                  $(FRAMEWORK_DIR)/kswr_cmp.cpp \
                  $(FRAMEWORK_DIR)/junit_reporter.cpp \
                  $(FRAMEWORK_DIR)/test_main.cpp
FRAMEWORK_OBJS := $(FRAMEWORK_SRCS:.cpp=.o)
FRAMEWORK_LIB  := $(FRAMEWORK_DIR)/libbwa_test_framework.a

FRAMEWORK_INCLUDES := $(INCLUDES) -I$(FRAMEWORK_DIR)

.PHONY: all framework unit integration clean libsais-objs
.SUFFIXES: .cpp .o

# Default: build everything.
all: framework unit integration $(EXE)

# -- Framework static lib ----------------------------------------------------

framework: $(FRAMEWORK_LIB)

$(FRAMEWORK_LIB): $(FRAMEWORK_OBJS)
	ar rcs $@ $^

$(FRAMEWORK_DIR)/%.o: $(FRAMEWORK_DIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(FRAMEWORK_INCLUDES) $< -o $@

# -- Doctest unit/integration binaries (populated in Task 12 + 13) ---------

unit: $(FRAMEWORK_LIB)
	@echo "unit target — no test sources yet (will be populated in Task 12)"

integration: $(FRAMEWORK_LIB)
	@echo "integration target — no test sources yet (will be populated in Task 13)"

# -- Legacy per-binary tests (MIGRATED ONE-BY-ONE IN FOLLOW-UP PRS) -------

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

# libsais object files live in the parent build tree and link into every
# test binary via $(LIBS) — they appear as bare paths inside the link
# command, NOT as make prerequisites, so make doesn't know to build them
# on demand. Delegate to the parent Makefile's LIBSAIS_OBJS rule
# (-DLIBSAIS_OPENMP / -O3 flag set stays in one place) and gate $(EXE)
# behind it via an order-only prerequisite — fires once per invocation
# without forcing the test binaries to relink. The parent's target paths
# are relative to its own dir, so we drop the leading `../` when
# invoking it.
libsais-objs:
	$(MAKE) -C .. ext/libsais/src/libsais.o ext/libsais/src/libsais64.o

$(EXE): | libsais-objs

fmi_test: fmi_test.o
	$(CXX) -o $@ $^ $(LIBS)

smem2_test: smem2_test.o
	$(CXX) -o $@ $^ $(LIBS)

sa2ref_test: sa2ref_test.o
	$(CXX) -o $@ $^ $(LIBS)

bwt_seed_strategy_test: bwt_seed_strategy_test.o
	$(CXX) -o $@ $^ $(LIBS)

xeonbsw: main_banded.o
	$(CXX) -o $@ $^ $(LIBS)

smem_lockstep_parity_test: smem_lockstep_parity_test.o
	$(CXX) -o $@ $^ $(LIBS)

header_insert_test:header_insert_test.o
	$(CXX) -o $@ $^ $(LIBS)

# Build test-local copies of packed_text.cpp / system.cpp into test/
# rather than ../src/. The parent build's CXXFLAGS (`-O3 -fpermissive`,
# `-DSAIS=1`, `-DLIBSAIS_OPENMP`) differ from this Makefile's
# (`-O2 -fopenmp -mtune=native -march=native`, `-DENABLE_PREFETCH`),
# and writing back into ../src/*.o would silently overwrite parent build
# artifacts with a different ABI/flag set — `ar rcs libbwa.a ...` would
# then archive the test-flagged objects on the next parent rebuild.
packed_text_test: packed_text_test.o packed_text_lib.o
	$(CXX) -o $@ $^ $(LIBS)

packed_text_test.o: ../src/packed_text.h
packed_text_lib.o: ../src/packed_text.cpp ../src/packed_text.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

system_test: system_test.o system_lib.o
	$(CXX) -o $@ $^ $(LIBS)

system_test.o: ../src/system.h
system_lib.o: ../src/system.cpp ../src/system.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@


clean:
	rm -f *.o $(FRAMEWORK_OBJS) $(FRAMEWORK_LIB) $(EXE) $(DOCTEST_BINS)

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
smem_lockstep_parity_test.o: ../src/FMI_search.h ../src/bntseq.h ../src/read_index_ele.h
smem_lockstep_parity_test.o: ../src/bwa.h ../src/bwt.h ../src/utils.h ../src/macro.h
header_insert_test.o: ../src/bwa.h
