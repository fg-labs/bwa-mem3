#!/usr/bin/env bash
# Negative controls for the manual gate: bash test/tools/disasm_identity_test.sh UNSTRIPPED_BINARY
# Uses a valid build for tool-failure tests, and synthetic symbols and listings for coverage,
# self-test and comparison failures.
set -euo pipefail
[ $# -eq 1 ] || {
    echo "usage: $0 UNSTRIPPED_BINARY" >&2
    exit 2
}
bin=$1
[[ $bin = /* ]] || bin=$PWD/$bin
gate=$(cd "$(dirname "$0")" && pwd)/disasm_identity.sh
tmp=$(mktemp -d)
trap 'chmod u+r "$tmp/unreadable" 2>/dev/null || true; rm -rf "$tmp"' EXIT
expect() {
    local want=$1 name=$2 rc=0
    shift 2
    "$@" > "$tmp/$name.log" 2>&1 || rc=$?
    if [ "$rc" != "$want" ]; then
        printf 'FAIL %s: expected exit %s, got %s\n' "$name" "$want" "$rc" >&2
        sed -n '1,25p' "$tmp/$name.log" >&2
        exit 1
    fi
    printf 'PASS %s (exit %s)\n' "$name" "$rc"
}
expect 0 self env DISASM_SELFTEST=1 bash "$gate" "$bin" "$bin"
expect 2 usage bash "$gate"
expect 2 nonexistent bash "$gate" "$tmp/missing" "$bin"
cp "$bin" "$tmp/unreadable"
chmod a-r "$tmp/unreadable"
if [ ! -r "$tmp/unreadable" ]; then
    expect 2 unreadable bash "$gate" "$tmp/unreadable" "$bin"
else
    echo 'SKIP unreadable (privileged user can still read it)'
fi
cp "$bin" "$tmp/stripped"
strip "$tmp/stripped" > "$tmp/strip.log" 2>&1
expect 2 stripped bash "$gate" "$tmp/stripped" "$bin"
printf 'not an executable\n' > "$tmp/text"
expect 2 nonbinary bash "$gate" "$tmp/text" "$bin"
mkdir "$tmp/tools"
real_nm=$(command -v nm)
real_objdump=$(command -v objdump)
# Tool failures, including the disassembly command failing after a successful version probe.
printf '#!/usr/bin/env bash\nexit 7\n' > "$tmp/tools/nm"
chmod +x "$tmp/tools/nm"
expect 2 nm-error env PATH="$tmp/tools:$PATH" bash "$gate" "$bin" "$bin"
ln -sf "$real_nm" "$tmp/tools/nm"
# The generated tool, not this test driver, must expand its positional argument.
# shellcheck disable=SC2016
printf '#!/usr/bin/env bash\nif [ "$1" = --version ]; then echo LLVM; exit 0; fi\nexit 7\n' > "$tmp/tools/objdump"
chmod +x "$tmp/tools/objdump"
expect 2 objdump-error env PATH="$tmp/tools:$PATH" bash "$gate" "$bin" "$bin"
printf '#!/usr/bin/env bash\nexit 7\n' > "$tmp/tools/objdump"
expect 2 objdump-version-error env PATH="$tmp/tools:$PATH" bash "$gate" "$bin" "$bin"
# shellcheck disable=SC2016 # Expansion belongs to the generated tool.
printf '#!/usr/bin/env bash\nif [ "$1" = --version ]; then echo LLVM; fi\nexit 0\n' > "$tmp/tools/objdump"
expect 2 empty-disassembly env PATH="$tmp/tools:$PATH" bash "$gate" "$bin" "$bin"
ln -sf "$real_objdump" "$tmp/tools/objdump"
rm "$tmp/tools/nm"
printf '#!/usr/bin/env bash\nexit 0\n' > "$tmp/tools/nm"
chmod +x "$tmp/tools/nm"
expect 2 all-symbols-absent env PATH="$tmp/tools:$PATH" bash "$gate" "$bin" "$bin"
# Exactly four defined stock symbols, including both required drivers. Return identical mock
# instructions for every function: ordinary comparison passes but DISASM_SELFTEST must exit 1.
rm "$tmp/tools/nm" "$tmp/tools/objdump"
printf '#!/usr/bin/env bash\n' > "$tmp/tools/nm"
for name in aln_search_batch bwa_cal_sa_reg_gap_group gap_bucket_grow gap_bucket_shrink; do
    printf "printf '00001000 T bwa_aln::%s()\\n'\n" "$name" >> "$tmp/tools/nm"
done
# shellcheck disable=SC2016 # Expansion belongs to the generated tool.
printf '#!/usr/bin/env bash\nif [ "$1" = --version ]; then echo LLVM; else echo "1000: ret"; fi\n' > "$tmp/tools/objdump"
chmod +x "$tmp/tools/nm" "$tmp/tools/objdump"
expect 0 minimal-core env PATH="$tmp/tools:$PATH" bash "$gate" "$bin" "$bin"
expect 1 selftest-failure env PATH="$tmp/tools:$PATH" DISASM_SELFTEST=1 bash "$gate" "$bin" "$bin"
# A lone unrelated helper must not make an empty/stripped comparison pass.
printf '#!/usr/bin/env bash\nprintf "00001000 T bwa_aln::gap_bucket_grow()\\n"\n' > "$tmp/tools/nm"
expect 2 insufficient-core env PATH="$tmp/tools:$PATH" bash "$gate" "$bin" "$bin"
# Comparisons of two different listings. The mock objdump prints the listing named by its binary
# (base or new) for the symbol it disassembles, so each side's instructions are chosen here.
cp "$bin" "$tmp/base"
cp "$bin" "$tmp/new"
printf '#!/usr/bin/env bash\nprintf "00001000 T bwa_aln::aln_search_batch()\\n00002000 T bwa_aln::bwa_cal_sa_reg_gap_group()\\n00003000 T bwa_aln::gap_bucket_grow()\\n00004000 T bwa_aln::gap_bucket_shrink()\\n"\n' > "$tmp/tools/nm"
# shellcheck disable=SC2016 # Expansion belongs to the generated tool.
printf '#!/usr/bin/env bash\n[ "$1" = --version ] && { echo LLVM; exit 0; }\nfor a; do case $a in --disassemble-symbols=*) s=${a#*=} ;; esac; done\nsed "s/SELF/$s/" "${!#}.listing"\n' > "$tmp/tools/objdump"
compare() { # expected exit, name, base listing, new listing (\n-separated instructions each)
    printf '%b' "$3" > "$tmp/base.listing"
    printf '%b' "$4" > "$tmp/new.listing"
    expect "$1" "$2" env PATH="$tmp/tools:$PATH" bash "$gate" "$tmp/base" "$tmp/new"
}
compare 1 base-new-mismatch '1000:\tret\n' '1000:\tnop\n'
# a branch keeps its offset into the function, whatever the function's address
compare 0 branch-moved '1000:\tb.ne\t0x1010 <SELF+0x10>\n' '2000:\tb.ne\t0x2010 <SELF+0x10>\n'
compare 1 branch-redirected '1000:\tb.ne\t0x1010 <SELF+0x10>\n' '1000:\tb.ne\t0x1014 <SELF+0x14>\n'
# a call keeps its callee: a template spelling of the same function is the same callee
compare 0 call-template-spelling '1000:\tbl\t0x5000 <_ZN7bwa_alnL11gap_job_hitEPNS_9gap_job_tEPKNS_5bwt_tE>\n' \
    '1000:\tbl\t0x6000 <_ZN7bwa_alnL11gap_job_hitINS_8AlnStockEEEiPNS_9gap_job_tEPKNS_5bwt_tE>\n'
compare 1 call-redirected '1000:\tcall\t5000 <malloc@plt>\n' '1000:\tcall\t5000 <calloc@plt>\n'
# an adrp page's address, its nearest-symbol annotation and its addend move with the layout
compare 0 adrp-addend '1000:\tadrp\tx0, 0x5000 <foo+0x10>\n1004:\tldr\tx1, [x0, #8]\n' \
    '1000:\tadrp\tx0, 0x6000 <bar+0x20>\n1004:\tldr\tx1, [x0, #16]\n'
# but a call can return anything in x0: an offset from it afterwards is compared, not normalised
compare 1 adrp-then-call '1000:\tadrp\tx0, 0x5000\n1004:\tbl\t0x5000 <foo>\n1008:\tldr\tx1, [x0, #8]\n' \
    '1000:\tadrp\tx0, 0x5000\n1004:\tbl\t0x5000 <foo>\n1008:\tldr\tx1, [x0, #16]\n'
