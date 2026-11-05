#!/usr/bin/env bash
# test/regression/input_fail_fast.sh
#
# Regression: `mem` reports an input it cannot read before it loads the index.
#
# The inputs are opened after the index load, which can take minutes on a large
# reference, so a typo in a read path used to cost a full load before the error.
# A stat() check now runs before the load; the real open and the reader's format
# sniff stay after it.
#
# Most legs give a MISSING INDEX on purpose: whichever check runs first decides
# the error, so "input error" means the input was rejected before the load and
# "index error" means it was left alone. Successful record comparisons and
# after-load permission checks use a tiny index; nothing touches the network.
#
#   1. A missing R1 or R2 (ENOENT/ENOTDIR) is reported (exit 1), not the index.
#      Under --meth too: the input is reported before the missing seed index.
#   2. A directory input is reported the same way.
#   3. The check does not open the input: with a FIFO that has no writer, the
#      index error must appear promptly (an open() would block forever; the run
#      is bounded by `timeout`).
#   4. Names kopen() handles itself are not checked -- http:// and ftp:// as
#      R1, http:// and '-' as R2, and under -p a missing second positional all
#      give the index error -- and
#      '-' (stdin), a `<cmd`, and a `<cmd` with leading whitespace still give
#      the plain-file records (which must be non-empty, so two header-only
#      outputs cannot pass).
#      Test-only stat interposition returns EOVERFLOW while real open/read
#      succeed: default/legacy readers preserve every R1/R2 record and defer
#      to a missing ordinary or methylation seed index.
#   5. An input that exists but is unreadable is not rejected before the load
#      (the index error appears); its own error stays after the load. Skipped
#      as root, which can read it anyway. Parent-directory EACCES is also
#      deferred, including with --meth; a missing R2 still precedes an
#      unreadable R1.
#   6. A missing input with -o exits 1 with the input error and creates no
#      output file (the -o file is only opened after the load). CI also runs
#      this script against its ASAN build, which fails on a memory error on
#      this early return. A leak of `opt` would not be caught: LeakSanitizer
#      still finds it reachable at exit.
#
# Inputs:
#   BWA_MEM3 — path to bwa-mem3 binary
#   FIXTURES — directory containing phix.fa and reads.fa (default: test/fixtures)
#   CC — C compiler (default: cc) for temporary Darwin/Linux interposition and
#        probe binaries. Missing CC or unsupported OS skips only injection;
#        compilation, probe or ineffective-injection errors fail the test.
# Leg 3 needs `timeout` (coreutils); without it the leg reports SKIP:.
set -euo pipefail

# Expected failures are captured explicitly below; diagnose unguarded setup errors.
trap 'status=$?; echo "FAIL: unexpected setup/check error at line $LINENO (exit $status)" >&2; exit "$status"' ERR

# Checked explicitly rather than with ${VAR:?}, which exits with a bare shell
# error (and before any ERR trap runs), so the failure would carry no FAIL: line.
if [[ -z ${BWA_MEM3:-} ]]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi
FIXTURES="${FIXTURES:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../fixtures" && pwd)}"

src_ref="$FIXTURES/phix.fa"
reads="$FIXTURES/reads.fa"
[[ -x "$BWA_MEM3" ]] || {
    echo "FAIL: binary not executable: $BWA_MEM3" >&2
    exit 1
}
[[ -s "$src_ref" ]] || {
    echo "FAIL: phix.fa missing: $src_ref" >&2
    exit 1
}
[[ -s "$reads" ]] || {
    echo "FAIL: reads.fa missing: $reads" >&2
    exit 1
}

mdir="$(mktemp -d)"
finished=0
cleanup() {
    local status=$?
    if [[ "$finished" -eq 0 ]]; then
        [[ "$status" -ne 0 ]] || status=1
        echo "FAIL: input fail-fast did not reach its final verdict (exit $status)" >&2
    fi
    chmod -R u+rwx "$mdir" 2> /dev/null || true
    rm -rf "$mdir" || true
    exit "$status"
}
trap cleanup EXIT
ref="$mdir/phix.fa"
cp "$src_ref" "$ref"
rc=0
"$BWA_MEM3" index "$ref" > "$mdir/index.log" 2>&1 || rc=$?
if [[ "$rc" -ne 0 ]]; then
    echo "FAIL: index phix.fa (exit $rc)"
    cat "$mdir/index.log"
    exit "$rc"
fi

fails=0
skips=0
fail() {
    echo "FAIL: $*"
    fails=$((fails + 1))
}
skip() {
    echo "SKIP: $*"
    skips=$((skips + 1))
}
show() { # $1 = rc, $2 = stderr file
    echo "      bwa-mem3 exit status $1; stderr tail:"
    tail -n 5 "$2" | sed 's/^/        /'
}

# Records of a SAM file without @ lines (the @PG CL: differs between runs).
records() { grep -v '^@' "$1" || true; }

noidx="$mdir/no-index/x.fa"
run() { # rest = mem args; sets rc, stderr in $mdir/e.err
    rc=0
    "$BWA_MEM3" mem "$@" > /dev/null 2> "$mdir/e.err" || rc=$?
}

# The input was rejected before the load: exit 1 and the input's message.
expect_input_error() { # $1 = label, $2 = expected message, rest = mem args
    local label="$1" want="$2"
    shift 2
    run "$@"
    if [[ "$rc" -eq 1 ]] && grep -qF "$want" "$mdir/e.err"; then
        echo "  ok: $label"
    else
        fail "$label: expected exit 1 with '$want'"
        show "$rc" "$mdir/e.err"
    fi
}
# The input was left alone and the (missing) index failed: exit 1, the index
# path named, no input error.
expect_index_error() { # $1 = label, rest = mem args (index is $noidx)
    local label="$1"
    shift
    run "$@"
    if [[ "$rc" -eq 1 ]] && grep -qF "$noidx" "$mdir/e.err" \
        && ! grep -q 'to open file `' "$mdir/e.err"; then
        echo "  ok: $label"
    else
        fail "$label: expected the index error (exit 1, naming $noidx), not an input error"
        show "$rc" "$mdir/e.err"
    fi
}

# --- 1. missing R1 / R2 ------------------------------------------------------
expect_input_error "missing R1 is reported before the index load" \
    "[E::main_mem] fail to open file \`$mdir/nope.1.fq'." "$noidx" "$mdir/nope.1.fq" "$reads"
expect_input_error "missing R2 is reported before the index load" \
    "[E::main_mem] failed to open file \`$mdir/nope.2.fq'." "$noidx" "$reads" "$mdir/nope.2.fq"
expect_input_error "--meth: missing R1 is reported before the missing seed index" \
    "[E::main_mem] fail to open file \`$mdir/nope.1.fq'." --meth "$noidx" "$mdir/nope.1.fq"
expect_input_error "--meth: missing R2 is reported before the missing seed index" \
    "[E::main_mem] failed to open file \`$mdir/nope.2.fq'." --meth "$noidx" "$reads" "$mdir/nope.2.fq"
expect_input_error "ENOTDIR in R1 is reported before the index load" \
    "[E::main_mem] fail to open file \`$reads/child'." "$noidx" "$reads/child"
expect_input_error "ENOTDIR in R2 is reported before the index load" \
    "[E::main_mem] failed to open file \`$reads/child'." "$noidx" "$reads" "$reads/child"
expect_index_error "--meth: with the inputs present, the seed index is reported" \
    --meth "$noidx" "$reads"

# --- 2. a directory ------------------------------------------------------------
mkdir "$mdir/a-dir"
expect_input_error "a directory input is reported before the index load" \
    "[E::main_mem] fail to open file \`$mdir/a-dir'." "$noidx" "$mdir/a-dir"
expect_input_error "a directory as R2 is reported before the index load, as R2" \
    "[E::main_mem] failed to open file \`$mdir/a-dir'." "$noidx" "$reads" "$mdir/a-dir"

# --- 3. the check does not open the input --------------------------------------
if command -v timeout > /dev/null 2>&1; then
    mkfifo "$mdir/nowriter.fifo"
    rc=0
    timeout 20 "$BWA_MEM3" mem "$noidx" "$mdir/nowriter.fifo" > /dev/null 2> "$mdir/e.err" || rc=$?
    if [[ "$rc" -eq 1 ]] && grep -qF "$noidx" "$mdir/e.err"; then
        echo "  ok: a FIFO with no writer is not opened before the load (index error at once)"
    elif [[ "$rc" -eq 124 ]]; then
        fail "a FIFO with no writer blocked the run for 20 s: the pre-load check opened it"
        show "$rc" "$mdir/e.err"
    else
        fail "a FIFO with no writer: expected the index error (exit 1)"
        show "$rc" "$mdir/e.err"
    fi
else
    skip "no 'timeout' to bound the FIFO run; the does-not-open check was not run"
fi

# --- 4. names the check leaves to kopen() ---------------------------------------
expect_index_error "http:// is not checked" "$noidx" "http://127.0.0.1:9/reads.fq"
expect_index_error "ftp:// is not checked" "$noidx" "ftp://127.0.0.1:9/reads.fq"
expect_index_error "R2 http:// is not checked" "$noidx" "$reads" "http://127.0.0.1:9/r2.fq"
expect_index_error "R2 '-' is not checked" "$noidx" "$reads" "-"
expect_index_error "R1 '-' is not checked" "$noidx" "-"
expect_index_error "a producer is not run before the index load" \
    "$noidx" "<touch '$mdir/producer-started'; cat '$reads'"
if [[ -e "$mdir/producer-started" ]]; then
    fail "preflight ran the producer before the index load"
fi
expect_index_error "-p does not check the ignored second positional" \
    -p "$noidx" "$reads" "$mdir/nope.2.fq"

rc=0
"$BWA_MEM3" mem "$ref" "$reads" > "$mdir/plain.sam" 2> "$mdir/plain.err" || rc=$?
if [[ "$rc" -ne 0 ]]; then
    fail "plain-file baseline run exited $rc"
    show "$rc" "$mdir/plain.err"
elif ! grep -qv '^@' "$mdir/plain.sam"; then
    fail "plain-file baseline produced no SAM records"
fi
records "$mdir/plain.sam" > "$mdir/plain.rec"
expected_reads=$(grep -c '^>' "$reads" || true)
plain_count=$(wc -l < "$mdir/plain.rec")
if [[ "$expected_reads" -eq 0 || "$plain_count" -ne "$expected_reads" ]]; then
    fail "plain-file baseline: expected $expected_reads records, got $plain_count"
else
    echo "  ok: plain-file baseline has all $plain_count records"
fi
special_input() { # $1 = label, $2 = the input argument (stdin is the reads)
    local label="$1" in="$2" rc=0
    "$BWA_MEM3" mem "$ref" "$in" < "$reads" > "$mdir/s.sam" 2> "$mdir/s.err" || rc=$?
    if [[ "$rc" -ne 0 ]]; then
        fail "$label: rejected or failed"
        show "$rc" "$mdir/s.err"
    elif ! grep -qv '^@' "$mdir/s.sam"; then
        fail "$label: produced no SAM records"
    elif ! diff -q "$mdir/plain.rec" <(records "$mdir/s.sam") > /dev/null; then
        fail "$label: output differs from the plain-file run"
    else
        echo "  ok: $label"
    fi
}
special_input "'-' (stdin) is accepted" "-"
special_input "a '<cmd' producer is accepted" "<cat '$reads'"
special_input "a '<cmd' with leading whitespace is accepted" "  <cat '$reads'"

# Metadata overflow need not prevent open/read (e.g. an unrepresentable inode).
# Keep the fault in a separate test library, with no production env-based seam.
# Reverting to rejection of every stat error must break the record comparisons
# and missing-index ordering checks below, for both input positions.
cp "$reads" "$mdir/stat-overflow.fa"
stat_src="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/input_stat_overflow.c"
stat_env=""
stat_flags=()
case "$(uname -s)" in
    Darwin)
        stat_flags=(-dynamiclib)
        stat_lib="$mdir/stat-overflow.dylib"
        stat_env="DYLD_INSERT_LIBRARIES=$mdir/stat-overflow.dylib"
        ;;
    Linux)
        stat_flags=(-shared -fPIC -ldl)
        stat_lib="$mdir/stat-overflow.so"
        stat_env="LD_PRELOAD=$mdir/stat-overflow.so"
        ;;
    *)
        stat_env=""
        skip "stat overflow injection needs Darwin or Linux"
        ;;
esac
if [[ -n "$stat_env" ]]; then
    if ! command -v "${CC:-cc}" > /dev/null 2>&1; then
        skip "stat overflow injection: C compiler unavailable: ${CC:-cc}"
        stat_env=""
    else
        rc=0
        "${CC:-cc}" "$stat_src" "${stat_flags[@]}" -o "$stat_lib" > "$mdir/compile.err" 2>&1 || rc=$?
        if [[ "$rc" -ne 0 ]]; then
            fail "stat overflow library compilation failed (exit $rc)"
            cat "$mdir/compile.err"
            stat_env=""
        else
            rc=0
            "${CC:-cc}" -DINPUT_STAT_PROBE "$stat_src" -o "$mdir/stat-probe" > "$mdir/compile.err" 2>&1 || rc=$?
            if [[ "$rc" -ne 0 ]]; then
                fail "stat overflow probe compilation failed (exit $rc)"
                cat "$mdir/compile.err"
                stat_env=""
            fi
        fi
    fi
fi
if [[ -n "$stat_env" ]]; then
    if env "$stat_env" "$mdir/stat-probe" "$mdir/stat-overflow.fa" > "$mdir/probe.fa" \
        && cmp -s "$reads" "$mdir/probe.fa"; then
        echo "  ok: stat EOVERFLOW while real open/read return every input byte"
    else
        fail "stat EOVERFLOW compatibility probe failed"
    fi
    # Distinct, matching-name mates derived from the reference expose accidental
    # R2-to-R1 substitution, which two copies of reads.fa cannot detect.
    sequence=$(grep -v '^>' "$src_ref" | tr -d '\r\n')
    : > "$mdir/paired-r1.fa"
    : > "$mdir/paired-r2.fa"
    for ((i = 0; i < 10; i++)); do
        start=$((i * 350))
        left=${sequence:start:100}
        right=${sequence:start+180:100}
        reverse=""
        for ((j = ${#right} - 1; j >= 0; j--)); do
            case "${right:j:1}" in
                A | a) reverse+=T ;;
                C | c) reverse+=G ;;
                G | g) reverse+=C ;;
                T | t) reverse+=A ;;
                *) reverse+=N ;;
            esac
        done
        [[ ${#left} -eq 100 && ${#reverse} -eq 100 && "$left" != "$reverse" ]]
        printf '>pair%s\n%s\n' "$i" "$left" >> "$mdir/paired-r1.fa"
        printf '>pair%s\n%s\n' "$i" "$reverse" >> "$mdir/paired-r2.fa"
    done
    for pos in SE R1 R2; do
        case "$pos" in
            SE)
                cp "$reads" "$mdir/stat-overflow.fa"
                inputs=("$mdir/stat-overflow.fa")
                ;;
            R1)
                cp "$mdir/paired-r1.fa" "$mdir/stat-overflow.fa"
                inputs=("$mdir/stat-overflow.fa" "$mdir/paired-r2.fa")
                ;;
            R2)
                cp "$mdir/paired-r2.fa" "$mdir/stat-overflow.fa"
                inputs=("$mdir/paired-r1.fa" "$mdir/stat-overflow.fa")
                ;;
        esac
        for reader in default legacy; do
            opts=()
            if [[ "$reader" == legacy ]]; then opts=(--legacy-reader); fi
            rc=0
            env "$stat_env" "$BWA_MEM3" mem ${opts[@]+"${opts[@]}"} "$ref" "${inputs[@]}" \
                > "$mdir/overflow.sam" 2> "$mdir/overflow.err" || rc=$?
            if ! grep -q '^test stat: EOVERFLOW$' "$mdir/overflow.err"; then
                fail "$pos/$reader: metadata interposition did not run"
            fi
            baseline_ok=1
            if [[ "$pos" == SE ]]; then
                expected_records="$mdir/plain.rec"
            else
                baseline_rc=0
                "$BWA_MEM3" mem ${opts[@]+"${opts[@]}"} "$ref" "$mdir/paired-r1.fa" "$mdir/paired-r2.fa" \
                    > "$mdir/paired.sam" 2> "$mdir/paired.err" || baseline_rc=$?
                if [[ "$baseline_rc" -ne 0 ]]; then
                    fail "$pos/$reader: paired baseline failed (exit $baseline_rc)"
                    show "$baseline_rc" "$mdir/paired.err"
                    baseline_ok=0
                fi
                records "$mdir/paired.sam" > "$mdir/paired.rec"
                paired_count=$(wc -l < "$mdir/paired.rec")
                if [[ "$paired_count" -ne 20 ]]; then
                    fail "paired baseline: expected 20 records, got $paired_count"
                    baseline_ok=0
                else
                    echo "  ok: paired baseline has all $paired_count records"
                fi
                expected_records="$mdir/paired.rec"
            fi
            if [[ "$baseline_ok" -ne 1 || "$rc" -ne 0 ]] || ! grep -qv '^@' "$mdir/overflow.sam" \
                || ! diff -q "$expected_records" <(records "$mdir/overflow.sam") > /dev/null; then
                fail "$pos/$reader: readable input rejected or records changed after stat EOVERFLOW"
                show "$rc" "$mdir/overflow.err"
            else
                echo "  ok: $pos/$reader: complete records unchanged after stat EOVERFLOW"
            fi
        done
        for mode in ordinary meth; do
            opts=()
            if [[ "$mode" == meth ]]; then opts=(--meth); fi
            rc=0
            env "$stat_env" "$BWA_MEM3" mem ${opts[@]+"${opts[@]}"} "$noidx" "${inputs[@]}" \
                > /dev/null 2> "$mdir/overflow.err" || rc=$?
            if [[ "$rc" -eq 1 ]] && grep -qF "$noidx" "$mdir/overflow.err" \
                && grep -q '^test stat: EOVERFLOW$' "$mdir/overflow.err" \
                && ! grep -q 'to open file `' "$mdir/overflow.err"; then
                echo "  ok: $pos/$mode: stat EOVERFLOW defers to the missing index"
            else
                fail "$pos/$mode: expected missing index, not stat EOVERFLOW input error"
                show "$rc" "$mdir/overflow.err"
            fi
        done
    done
fi

# --- 5. an existing but unreadable input is left to the open after the load ----
if [[ "$(id -u)" -eq 0 ]]; then
    skip "running as root, which can read a mode-000 file; unreadable-input check not run"
else
    cp "$reads" "$mdir/unreadable.fa"
    chmod 000 "$mdir/unreadable.fa"
    expect_index_error "an unreadable input is not rejected before the load" \
        "$noidx" "$mdir/unreadable.fa"
    expect_input_error "an unreadable input is still reported after the load" \
        "[E::main_mem] fail to open file \`$mdir/unreadable.fa'." "$ref" "$mdir/unreadable.fa"
    expect_index_error "an unreadable R2 is not rejected before the load" \
        "$noidx" "$reads" "$mdir/unreadable.fa"
    expect_input_error "an unreadable R2 is still reported after the load" \
        "[E::main_mem] failed to open file \`$mdir/unreadable.fa'." "$ref" "$reads" "$mdir/unreadable.fa"
    expect_input_error "missing R2 precedes the unreadable R1 open" \
        "[E::main_mem] failed to open file \`$mdir/nope.2.fq'." "$noidx" "$mdir/unreadable.fa" "$mdir/nope.2.fq"
    mkdir "$mdir/unsearchable"
    cp "$reads" "$mdir/unsearchable/reads.fa"
    chmod 000 "$mdir/unsearchable"
    expect_index_error "parent EACCES in R1 is left to the open after the load" \
        "$noidx" "$mdir/unsearchable/reads.fa"
    expect_index_error "parent EACCES in R2 is left to the open after the load" \
        "$noidx" "$reads" "$mdir/unsearchable/reads.fa"
    expect_index_error "--meth: parent EACCES is left until after the seed index" \
        --meth "$noidx" "$mdir/unsearchable/reads.fa"
    expect_input_error "parent EACCES in R1 is still reported after the load" \
        "[E::main_mem] fail to open file \`$mdir/unsearchable/reads.fa'." "$ref" "$mdir/unsearchable/reads.fa"
    expect_input_error "parent EACCES in R2 is still reported after the load" \
        "[E::main_mem] failed to open file \`$mdir/unsearchable/reads.fa'." "$ref" "$reads" "$mdir/unsearchable/reads.fa"
    chmod 700 "$mdir/unsearchable"
fi

# --- 6. a missing input with -o ------------------------------------------------
expect_input_error "a missing input with -o exits cleanly" \
    "[E::main_mem] fail to open file \`$mdir/nope.1.fq'." \
    -o "$mdir/out.sam" "$noidx" "$mdir/nope.1.fq"
if [[ -e "$mdir/out.sam" ]]; then
    fail "a missing input with -o created the output file"
fi

finished=1
if [[ "$fails" -ne 0 ]]; then
    echo "FAIL: input fail-fast ($fails failure(s), $skips skipped)"
    exit 1
fi
if [[ "$skips" -ne 0 ]]; then
    echo "PASS: input fail-fast ($skips check(s) skipped; see SKIP: above)"
else
    echo "PASS: a missing input or a directory is reported before the index load; FIFOs, stdin, <cmd, URLs, -p and unreadable files are left to the open after it"
fi
