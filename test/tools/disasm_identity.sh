#!/usr/bin/env bash
# test/tools/disasm_identity.sh BASELINE_BINARY NEW_BINARY
#
# Manual pre-release check (CI runs only its negative controls, disasm_identity_test.sh): build
# both revisions with the same compiler and flags. Compare stock aln instructions, normalising only
# addresses, PC-relative addresses and relocation addends. Branch and call targets are compared by
# symbol (a branch by its offset into the function), and memory-operand displacements and other
# immediates must stay identical.
#
# Positive control: -DALN_DISASM_POSITIVE_CONTROL adds 16 bytes to gap_slot_core_t and must fail:
#   make EXTRA_CXXFLAGS=-DALN_DISASM_POSITIVE_CONTROL ... ; disasm_identity.sh BASE NEW -> exit 1
#
# Exit 0: identical. Exit 1: instructions differ or the self-test fails. Exit 2: invalid inputs,
# tool errors, missing/ambiguous symbols, empty disassembly or insufficient coverage. Stripped
# binaries are not supported. DISASM_SELFTEST=1 also compares the batch and group drivers with
# each other; the self-test must detect that they differ.
set -euo pipefail
fail() {
    echo "ERROR: $*" >&2
    exit 2
}
[ $# -eq 2 ] || fail "usage: $0 BASELINE NEW"
base=$1 new=$2
# Do not let a relative filename beginning with '-' become a tool option.
[[ $base = /* ]] || base=$PWD/$base
[[ $new = /* ]] || new=$PWD/$new
tmp=$(mktemp -d) || fail "cannot create temporary directory"
trap 'rm -rf "$tmp"' EXIT
for side in base new; do
    bin=${!side}
    [ -f "$bin" ] && [ -r "$bin" ] || fail "$bin is not a readable file"
    nm "$bin" > "$tmp/$side.nm" 2> "$tmp/nm.err" || fail "nm failed on $bin: $(< "$tmp/nm.err")"
    nm -C "$bin" > "$tmp/$side.demangled" 2> "$tmp/nm.err" || fail "nm -C failed on $bin: $(< "$tmp/nm.err")"
done
objdump --version > "$tmp/objdump.version" 2>&1 || fail "objdump --version failed"
llvm=0
if grep -qi llvm "$tmp/objdump.version"; then llvm=1; fi
# LTO can inline the exported job/slot functions into the group driver, and can remove the
# unused reference search. Require BOTH drivers (search/width loop and batch/dedup) plus at
# least two other stock functions, normally the cold hit/stack helpers that survive LTO.
MIN_COMPARED=4
BATCH='bwa_aln::aln_search_batch('
GROUP='bwa_aln::bwa_cal_sa_reg_gap_group('
# Old|templated demangled prefixes. Either spelling is accepted on either side so that the gate
# also works for two post-template builds (including a positive-control build).
FUNCS=(
    'bwa_aln::gap_slots_run('
    'bwa_aln::gap_job_start('
    'bwa_aln::gap_job_run('
    'bwa_aln::bwt_match_gap('
    "$BATCH"
    'bwa_aln::search_reads(|bwa_aln::search_reads<bwa_aln::AlnStock>('
    'bwa_aln::gap_job_hit(|bwa_aln::gap_job_hit<bwa_aln::AlnStock>('
    'bwa_aln::gap_bucket_grow('
    'bwa_aln::gap_bucket_shrink(|bwa_aln::gap_bucket_shrink<bwa_aln::gap_entry_t>('
    "$GROUP|bwa_aln::bwa_cal_sa_reg_gap_group<bwa_aln::AlnStock>("
)
symbol() { # side, old prefix, templated prefix -> address (only defined text, no lambdas/cold clones)
    mawk -v a="$2" -v b="$3" '
        $2 ~ /^[TtWw]$/ && (index($0, a) || index($0, b)) &&
        $0 !~ /::\$_|clERN|lambda|\.cold|\[clone/ { print $1 }
    ' "$tmp/$1.demangled"
}
dis() { # side, symbol address -> normalised listing on stdout
    local side=$1 addr=$2 mang
    mang=$(mawk -v a="$addr" '$1 == a && $2 ~ /^[TtWw]$/ {print $3; exit}' "$tmp/$side.nm")
    [ -n "$mang" ] || fail "cannot resolve symbol at $addr in $side"
    if [ "$llvm" = 1 ]; then
        objdump -d --no-show-raw-insn --disassemble-symbols="$mang" "${!side}" > "$tmp/raw" 2> "$tmp/objdump.err" || fail "objdump failed on $mang: $(< "$tmp/objdump.err")"
    else
        objdump -d --no-show-raw-insn --disassemble="$mang" "${!side}" > "$tmp/raw" 2> "$tmp/objdump.err" || fail "objdump failed on $mang: $(< "$tmp/objdump.err")"
    fi
    # Branch and call targets are compared by symbol, not address: objdump's <symbol+offset>
    # annotation is kept, with this function's own name as SELF (a branch keeps its offset into
    # the function) and any other symbol reduced to its demangled name without template
    # arguments, parameters or return type (a call keeps its callee, under either spelling).
    grep -oE '<[^<>]+>' "$tmp/raw" | sed -E -e 's/^<//' -e 's/>$//' -e 's/\+0x[0-9a-f]+$//' | sort -u > "$tmp/names" || true
    c++filt < "$tmp/names" > "$tmp/names.demangled" 2> "$tmp/cxxfilt.err" || fail "c++filt failed: $(< "$tmp/cxxfilt.err")"
    paste "$tmp/names" "$tmp/names.demangled" > "$tmp/names.map"
    # Only instruction lines, not file/section/symbol headings. An empty result is an error,
    # including objdump implementations that warn about a missing symbol but return success.
    sed -nE '/^[[:space:]]*[0-9a-f]+:[[:space:]]+/p' "$tmp/raw" \
        | sed -E -e 's/^[[:space:]]*[0-9a-f]+:[[:space:]]*//' -e 's/[[:space:]]*;.*$//' -e 's/[[:space:]]+#([[:space:]]|$).*$//' \
        | mawk -v self="$mang" -v mapf="$tmp/names.map" '
            function canon(s,   out, i, c, d, sfx, n, t) { # demangled name -> scope::name
                sfx = ""
                if ((i = index(s, " [clone ")) > 0) { sfx = substr(s, i); s = substr(s, 1, i - 1) }
                out = ""; d = 0
                for (i = 1; i <= length(s); i++) {
                    c = substr(s, i, 1)
                    if (c == "<" || c == "(") d++
                    else if (c == ">" || c == ")") { if (d > 0) d-- }
                    else if (d == 0) out = out c
                }
                n = split(out, t, " ")
                while (n > 1 && t[n] ~ /^(const|volatile|&|&&|noexcept)$/) n-- # drop a method qualifier
                return t[n] sfx
            }
            BEGIN { while ((getline line < mapf) > 0) { split(line, kv, "\t"); dem[kv[1]] = kv[2] } }
            { m = $1
                # PC-relative targets (branches, calls, adr/adrp pages) move with unrelated code.
                tgt = ""
                if (m ~ /^(b|bl|br|blr|cbz|cbnz|tbz|tbnz|adr|adrp|call|callq|jmp|jmpq)$/ || m ~ /^b\./ || m ~ /^j[a-z]+$/) {
                  # (an adr/adrp page is annotated with whatever symbol lies nearest: layout, not a target)
                  if (m !~ /^adrp?$/ && match($0, /<[^<>]+>/)) {
                    tgt = substr($0, RSTART + 1, RLENGTH - 2); off = ""
                    if (match(tgt, /\+0x[0-9a-f]+$/)) { off = substr(tgt, RSTART); tgt = substr(tgt, 1, RSTART - 1) }
                    if (tgt == self) tgt = "<SELF" off ">"
                    else tgt = "<" canon((tgt in dem) ? dem[tgt] : tgt) ">"
                  }
                  gsub(/(0x)?[0-9a-f][0-9a-f][0-9a-f][0-9a-f]+/, "ADDR")
                }
                gsub(/<[^<>]*>/, "")
                gsub(/0x[0-9a-f]+\(%rip\)/, "ADDR(%rip)"); gsub(/rip \+ 0x[0-9a-f]+/, "rip + ADDR")
                # arm64 low 12-bit relocation addends for registers still holding an adrp page.
                ops = $0; sub(/^[^ \t]+[ \t]+/, "", ops); split(ops, a, ",")
                if (m == "adrp") { pg[a[1]] = 1; print; next }
                # A call may return anything in the caller-saved registers (AAPCS64 x0-x18, and
                # x30 holds the return address): none of them still holds its page afterwards.
                if (m == "bl" || m == "blr") { for (i = 0; i <= 18; i++) delete pg["x" i]; delete pg["x30"] }
                if (m == "add" && (a[2] ~ /^ *x[0-9]+$/)) { r = a[2]; gsub(/ /, "", r); if (r in pg) sub(/#0x[0-9a-f]+|#[0-9]+/, "#LO12") }
                if (match($0, /\[x[0-9]+, #/)) { r = substr($0, RSTART + 1); sub(/,.*/, "", r); if (r in pg) sub(/\[x[0-9]+, #(0x[0-9a-f]+|[0-9]+)/, "[" r ", #LO12") }
                d = a[1]; gsub(/ /, "", d)
                if ((d in pg) && m !~ /^(st|cmp|cmn|tst|ccmp)/) delete pg[d]
                print $0 tgt }'
}
rc=0
n_compared=0
batch_seen=0
group_seen=0
for spec in "${FUNCS[@]}"; do
    f=${spec%%|*}
    g=${spec#*|}
    a=$(symbol base "$f" "$g") || fail "cannot read baseline symbols"
    b=$(symbol new "$f" "$g") || fail "cannot read new symbols"
    if [ -z "$a" ] && [ -z "$b" ]; then
        echo "ABSENT $f (inlined or unused in both builds)"
        continue
    fi
    [ -n "$a" ] && [ -n "$b" ] || fail "missing $f on one side"
    [[ $a != *$'\n'* && $b != *$'\n'* ]] || fail "ambiguous symbol $f"
    dis base "$a" > "$tmp/a" || fail "cannot disassemble $f in baseline"
    dis new "$b" > "$tmp/b" || fail "cannot disassemble $f in new binary"
    [ -s "$tmp/a" ] && [ -s "$tmp/b" ] || fail "empty disassembly of $f"
    n_compared=$((n_compared + 1))
    if [ "$f" = "$BATCH" ]; then
        batch_seen=1
        cp "$tmp/a" "$tmp/batch"
    fi
    if [ "$f" = "$GROUP" ]; then
        group_seen=1
        cp "$tmp/a" "$tmp/group"
    fi
    if cmp -s "$tmp/a" "$tmp/b"; then
        echo "IDENTICAL $f ($(wc -l < "$tmp/a") instructions)"
    else
        echo "DIFFERS $f"
        diff "$tmp/a" "$tmp/b" > "$tmp/d" || true
        sed -n 1,20p "$tmp/d"
        rc=1
    fi
done
if [ "$n_compared" -lt "$MIN_COMPARED" ] || [ "$batch_seen" = 0 ] || [ "$group_seen" = 0 ]; then
    fail "only $n_compared functions compared (need $MIN_COMPARED, including batch and group drivers); stripped or unsupported symbol layout?"
fi
echo "compared $n_compared functions"
if [ "${DISASM_SELFTEST:-0}" = 1 ]; then
    if cmp -s "$tmp/batch" "$tmp/group"; then
        echo "SELFTEST FAIL: the comparison cannot tell two functions apart"
        rc=1
    else echo "SELFTEST ok: different functions compare different"; fi
fi
exit "$rc"
