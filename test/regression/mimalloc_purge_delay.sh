#!/usr/bin/env bash
# test/regression/mimalloc_purge_delay.sh
#
# Regression: `bwa-mem3 mem` turns off mimalloc's page purging by default
# (purge_delay = -1), an explicit MIMALLOC_PURGE_DELAY still wins, a value
# mimalloc cannot parse is reported rather than silently dropped, and the
# setting never changes the alignments.
#
# mimalloc v3 decommits freed pages after `purge_delay` ms (1000 by default).
# `mem` frees and reallocates the same large per-batch buffers every batch, so
# that purge turns buffer reuse into repeated page faults and system time;
# `mem` therefore defaults it to -1 (never purge) through
# mi_option_set_default, which yields to a value mimalloc already read from
# the environment. The `[M::main] mimalloc purge delay: <N> ms` stderr line
# reports the value in effect (read back with mi_option_get after the default
# is applied), which is what the checks below observe.
#
# Checks:
#   1. default `mem` runs with purge delay -1
#   2. MIMALLOC_PURGE_DELAY=250 overrides it to 250, and so do the legacy
#      alias MIMALLOC_RESET_DELAY=250 and the lower-case spelling
#      mimalloc_purge_delay=250 (bwa-mem3 reads the same names for its
#      warning, so the two must agree)
#   3. an unparseable value (`1s`) falls back to -1 WITH a `[W::main]`
#      warning naming it, including under a lower-case variable name (mimalloc
#      matches names case-insensitively), while an explicit -1 and a boolean
#      word mimalloc does accept (`off` -> 0) produce no warning
#   4. the SAM is byte-identical between the default (never purge) and
#      MIMALLOC_PURGE_DELAY=0 (purge immediately, the opposite extreme), over
#      several batches (`-K`), so buffers really are freed and reused between
#      batches in both arms
#   5. `index` does not print the purge-delay line. This shows the `mem`
#      reporting path is not reached from `index`; it cannot show that the
#      mi_option_set_default call itself stays in the `mem` branch, since a
#      moved call would print nothing either.
#
# A binary built with USE_MIMALLOC=0, or one whose mimalloc does not override
# malloc, has no purge delay to set: the script reports SKIP, unless
# PURGE_MIMALLOC_REQUIRED is set, in which case it FAILs. CI sets it for the
# canonical build, which must run on mimalloc, so losing the malloc override
# there cannot pass as a skip.
#
# Fixture (deterministic, no PRNG): single-end reads sliced from the committed
# phix.fa.
#
# Inputs (env vars):
#   BWA_MEM3         — path to the bwa-mem3 binary under test
#   PURGE_PHIX_FA    — path to test/fixtures/phix.fa (the reference source)
#   PURGE_MIMALLOC_REQUIRED — optional; when non-empty, a binary that does not
#                      run on mimalloc is a FAIL rather than a SKIP
#   PURGE_WORK_DIR   — parent directory for fixture-private intermediates; each
#                      invocation works in its own subdirectory, removed on exit

set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${PURGE_PHIX_FA:?PURGE_PHIX_FA must be set}"
: "${PURGE_WORK_DIR:?PURGE_WORK_DIR must be set}"

# The checks below are about the defaults, so a caller's own setting must not
# leak into the "default" runs. mimalloc matches these names case-insensitively,
# so clear both spellings.
unset MIMALLOC_PURGE_DELAY MIMALLOC_RESET_DELAY mimalloc_purge_delay mimalloc_reset_delay

# Captured before grepping: piping `version` into `grep -q` would turn a crashing
# binary (or a SIGPIPE under pipefail) into a SKIP rather than a FAIL.
if ! version_out="$("$BWA_MEM3" version 2>&1)"; then
    echo "FAIL: $BWA_MEM3 version exited nonzero" >&2
    printf '%s\n' "$version_out" >&2
    exit 1
fi
if ! grep -Eq '^mimalloc [0-9]+\.[0-9]+\.[0-9]+ \(active\)$' <<< "$version_out"; then
    if [[ -n ${PURGE_MIMALLOC_REQUIRED:-} ]]; then
        echo "FAIL: PURGE_MIMALLOC_REQUIRED is set but $BWA_MEM3 does not run on mimalloc" >&2
        printf '%s\n' "$version_out" >&2
        exit 1
    fi
    echo "SKIP: $BWA_MEM3 does not run on mimalloc; no purge delay to check"
    exit 0
fi

# A private subdirectory per invocation, so concurrent runs sharing
# PURGE_WORK_DIR never read each other's reads, logs or SAMs.
mkdir -p "$PURGE_WORK_DIR"
PURGE_WORK_DIR="$(mktemp -d "$PURGE_WORK_DIR/run.XXXXXX")"
trap 'rm -rf "$PURGE_WORK_DIR"' EXIT

ref="$PURGE_WORK_DIR/phix.fa"
cp "$PURGE_PHIX_FA" "$ref"

# --- Build a single-end FASTQ by slicing phix (deterministic, no PRNG). ---
seq="$(grep -v '^>' "$ref" | tr -d '\n' | tr 'acgt' 'ACGT')"
fq="$PURGE_WORK_DIR/reads.fq"
: > "$fq"
L=100
q="$(printf 'I%.0s' $(seq 1 "$L"))"
n_reads=0
for off in $(seq 0 150 5200); do
    s="${seq:$off:$L}"
    [ "${#s}" -eq "$L" ] || continue
    printf '@read%d\n%s\n+\n%s\n' "$n_reads" "$s" "$q" >> "$fq"
    n_reads=$((n_reads + 1))
done
if [ "$n_reads" -lt 10 ]; then
    echo "FAIL: fixture built only $n_reads reads (phix slicing regression?)" >&2
    exit 1
fi

# --- Check 5: `index` does not reach the `mem` reporting path. ---
index_log="$PURGE_WORK_DIR/index.log"
if ! "$BWA_MEM3" index "$ref" > "$index_log" 2>&1; then
    echo "FAIL: bwa-mem3 index exited nonzero" >&2
    cat "$index_log" >&2
    exit 1
fi
if grep -q 'mimalloc purge delay' "$index_log"; then
    echo "FAIL: \`index\` reports a mimalloc purge delay; the default is meant for \`mem\` only" >&2
    cat "$index_log" >&2
    exit 1
fi

# run_mem NAME [VAR=VALUE ...] [-- MEM_ARG ...]
# Run `mem` with the given environment assignments and extra `mem` arguments,
# writing the SAM to NAME.sam and stderr to NAME.log. Exits on a nonzero status.
# The `--` keeps the two lists apart: assignments go to `env`, arguments to `mem`.
run_mem() {
    local name="$1"
    shift
    local -a env_args=() mem_args=()
    while (($# > 0)) && [ "$1" != "--" ]; do
        env_args+=("$1")
        shift
    done
    if (($# > 0)); then
        shift
        mem_args=("$@")
    fi
    if ! env ${env_args[@]+"${env_args[@]}"} "$BWA_MEM3" mem ${mem_args[@]+"${mem_args[@]}"} "$ref" "$fq" \
        > "$PURGE_WORK_DIR/$name.sam" 2> "$PURGE_WORK_DIR/$name.log"; then
        echo "FAIL: bwa-mem3 mem ($name) exited nonzero" >&2
        cat "$PURGE_WORK_DIR/$name.log" >&2
        exit 1
    fi
}

# Assert that run $1 reported purge delay $2.
expect_delay() {
    local name="$1" want="$2" got
    got="$(sed -nE 's/^\[M::main\] mimalloc purge delay: (-?[0-9]+) ms.*/\1/p' "$PURGE_WORK_DIR/$name.log")"
    if [ "$got" != "$want" ]; then
        echo "FAIL: $name: expected purge delay '$want', got '${got:-<no line>}'" >&2
        cat "$PURGE_WORK_DIR/$name.log" >&2
        exit 1
    fi
    echo "ok: $name runs with purge delay $want"
}

# Assert that run $1 did (want=yes) or did not (want=no) warn about an ignored
# purge-delay variable.
expect_warning() {
    local name="$1" want="$2" have=no
    if grep -qi '^\[W::main\] ignoring MIMALLOC_' "$PURGE_WORK_DIR/$name.log"; then
        have=yes
    fi
    if [ "$have" != "$want" ]; then
        echo "FAIL: $name: expected ignored-value warning: $want, got: $have" >&2
        cat "$PURGE_WORK_DIR/$name.log" >&2
        exit 1
    fi
}

# --- Check 1: the default. ---
run_mem default
expect_delay default -1
expect_warning default no

# --- Check 2: both environment names override it. ---
run_mem env_override MIMALLOC_PURGE_DELAY=250
expect_delay env_override 250
expect_warning env_override no
run_mem legacy_override MIMALLOC_RESET_DELAY=250
expect_delay legacy_override 250
expect_warning legacy_override no
run_mem lowercase_override mimalloc_purge_delay=250
expect_delay lowercase_override 250
expect_warning lowercase_override no

# --- Check 3: values mimalloc cannot parse are reported, not silently dropped. ---
run_mem malformed MIMALLOC_PURGE_DELAY=1s
expect_delay malformed -1
expect_warning malformed yes
if ! grep -q "MIMALLOC_PURGE_DELAY='1s'" "$PURGE_WORK_DIR/malformed.log"; then
    echo "FAIL: malformed: the warning does not name the ignored value" >&2
    cat "$PURGE_WORK_DIR/malformed.log" >&2
    exit 1
fi
run_mem legacy_malformed MIMALLOC_RESET_DELAY=1s
expect_delay legacy_malformed -1
expect_warning legacy_malformed yes
# mimalloc matches the name case-insensitively, so it reads (and drops) this too.
run_mem lowercase_malformed mimalloc_purge_delay=1s
expect_delay lowercase_malformed -1
expect_warning lowercase_malformed yes
if ! grep -q "mimalloc_purge_delay='1s'" "$PURGE_WORK_DIR/lowercase_malformed.log"; then
    echo "FAIL: lowercase_malformed: the warning does not name the variable as spelled" >&2
    cat "$PURGE_WORK_DIR/lowercase_malformed.log" >&2
    exit 1
fi
run_mem explicit_never MIMALLOC_PURGE_DELAY=-1
expect_delay explicit_never -1
expect_warning explicit_never no
run_mem boolean_word MIMALLOC_PURGE_DELAY=off
expect_delay boolean_word 0
expect_warning boolean_word no

# --- Check 4: the purge setting never changes the alignments. ---
# A small -K splits the fixture into several batches, so each arm frees the
# per-batch buffers and allocates them again; one batch would never exercise
# the reuse this default exists for.
batch_args=(-K 1000)
run_mem batched_never -- "${batch_args[@]}"
expect_delay batched_never -1
run_mem batched_purge_now MIMALLOC_PURGE_DELAY=0 -- "${batch_args[@]}"
expect_delay batched_purge_now 0
n_batches="$(grep -c 'Processed' "$PURGE_WORK_DIR/batched_never.log" || true)"
if [ "$n_batches" -lt 2 ]; then
    echo "FAIL: -K ${batch_args[1]} gave $n_batches batch(es); the identity check needs several" >&2
    cat "$PURGE_WORK_DIR/batched_never.log" >&2
    exit 1
fi
if [ "$(grep -vc '^@' "$PURGE_WORK_DIR/batched_never.sam")" -ne "$n_reads" ]; then
    echo "FAIL: batched run emitted the wrong number of records for $n_reads reads" >&2
    exit 1
fi
if ! diff -u "$PURGE_WORK_DIR/batched_never.sam" "$PURGE_WORK_DIR/batched_purge_now.sam" > "$PURGE_WORK_DIR/purge.diff"; then
    echo "FAIL: SAM differs between purge delay -1 and 0:" >&2
    head -40 "$PURGE_WORK_DIR/purge.diff" >&2
    exit 1
fi
echo "ok: SAM byte-identical between purge delay -1 and 0 over $n_batches batches"

echo "PASS: mem defaults mimalloc purge delay to -1, env overrides it, unparseable values warn, and output is byte-identical ($n_reads reads, phix)"
