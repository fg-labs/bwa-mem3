#!/usr/bin/env bash
# Negative controls for the manual gate: bash test/tools/disasm_identity_test.sh UNSTRIPPED_BINARY
# Uses a valid build for tool-failure tests, and synthetic symbols for coverage/self-test failures.
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
