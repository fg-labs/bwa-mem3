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
# "index error" means it was left alone. Nothing here loads a real index except
# the stdin / `<cmd` legs, and nothing touches the network.
#
#   1. A missing R1 or R2 is reported (usual message, exit 1), not the index.
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
#   5. An input that exists but is unreadable is not rejected before the load
#      (the index error appears); its own error stays after the load. Skipped
#      as root, which can read it anyway.
#   6. A missing input with -o exits 1 with the input error and creates no
#      output file (the -o file is only opened after the load). CI also runs
#      this script against its ASAN build, which fails on a memory error on
#      this early return. A leak of `opt` would not be caught: LeakSanitizer
#      still finds it reachable at exit.
#
# Inputs:
#   BWA_MEM3 — path to bwa-mem3 binary
#   FIXTURES — directory containing phix.fa and reads.fa (default: test/fixtures)
# Leg 3 needs `timeout` (coreutils); without it the leg reports SKIP:.
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
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
trap 'chmod -R u+rwx "$mdir" 2> /dev/null || true; rm -rf "$mdir"' EXIT
ref="$mdir/phix.fa"
cp "$src_ref" "$ref"
if ! "$BWA_MEM3" index "$ref" > /dev/null 2>&1; then
    echo "FAIL: index phix.fa"
    exit 1
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
fi

# --- 6. a missing input with -o ------------------------------------------------
expect_input_error "a missing input with -o exits cleanly" \
    "[E::main_mem] fail to open file \`$mdir/nope.1.fq'." \
    -o "$mdir/out.sam" "$noidx" "$mdir/nope.1.fq"
if [[ -e "$mdir/out.sam" ]]; then
    fail "a missing input with -o created the output file"
fi

if [[ "$fails" -ne 0 ]]; then
    echo "FAIL: input fail-fast ($fails failure(s), $skips skipped)"
    exit 1
fi
if [[ "$skips" -ne 0 ]]; then
    echo "PASS: input fail-fast ($skips check(s) skipped; see SKIP: above)"
else
    echo "PASS: a missing input or a directory is reported before the index load; FIFOs, stdin, <cmd, URLs, -p and unreadable files are left to the open after it"
fi
