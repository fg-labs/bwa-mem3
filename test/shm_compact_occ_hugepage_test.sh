#!/usr/bin/env bash
# Huge-page legs of shm_compact_occ_test.sh, registered separately so an
# unavailable hugetlb stage does not hide passing normal-layout coverage.
# The incompatible-options refusal must still run on every platform.

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR/.."

BIN=${BIN:-./bwa-mem3}
PREFIX=test/fixtures/phix.fa

if [[ ! -x "$BIN" ]]; then
    echo "FAIL: $BIN not built. Run 'make -j4' first." >&2
    exit 2
fi

# shellcheck source=test/framework/shm_registry_guard.sh
source "$SCRIPT_DIR/framework/shm_registry_guard.sh"
# shellcheck source=test/framework/shm_test_helpers.sh
source "$SCRIPT_DIR/framework/shm_test_helpers.sh"
shm_registry_guard "$BIN" shm_compact_occ_hugepage_test phix.fa

need_index=0
for ext in .amb .ann .bwt.2bit.64 .pac; do
    if [[ ! -s "${PREFIX}${ext}" ]]; then
        need_index=1
        break
    fi
done
if [[ "$need_index" -eq 1 ]]; then
    echo "[setup] Building phix index..."
    "$BIN" index "$PREFIX" > /dev/null || {
        echo "FAIL: index failed for '$PREFIX'" >&2
        exit 1
    }
fi

WORK=$(mktemp -d -t shm_compact_hp.XXXXXX)
trap 'shm_test_exit_cleanup "$?"' EXIT
unset BWA3_COMPACT_OCC
shm_test_destroy "setup"

fail() {
    echo "FAIL: $1" >&2
    shift
    for f in "$@"; do
        echo "----- $f (last 20 lines) -----" >&2
        tail -20 "$f" >&2
    done
    exit 1
}
list_registry() {
    "$BIN" shm -l > "$WORK/list.out" 2> "$WORK/list.err" \
        || fail "shm -l failed" "$WORK/list.err" "$WORK/list.out"
    # Global destroy cannot select this stage alone. On seeing a foreign row,
    # preserve every stage, including on EXIT; remove only our temporary files.
    if ! awk -F'\t' '$1 != "phix.fa" { exit 1 }' "$WORK/list.out"; then
        trap 'shm_test_temp_cleanup "$?"' EXIT
        fail "G: refusing to clear an unrelated registry entry" "$WORK/list.out"
    fi
}
require_default_layout() {
    list_registry
    awk -F'\t' '
        $1 == "phix.fa" {
            count++
            if (NF == 2 && $2 ~ /^[0-9]+$/ && $2 > 0) valid++
        }
        END { exit !(count == 1 && valid == 1) }
    ' "$WORK/list.out" || fail "shm -l: expected exactly one default row for phix.fa" "$WORK/list.out"
}
require_empty_registry() {
    list_registry
    [[ ! -s "$WORK/list.out" ]] || fail "shm -l: expected empty registry" "$WORK/list.out"
}

echo "[G] shm --huge-pages --compact-occ is refused"
# A 1 GB stage holds the default checkpoint layout only. The refusal precedes
# every host check (Linux, AMD, the hugetlb pool), so test it before trying one.
refusal_status=0
"$BIN" shm --huge-pages --compact-occ "$PREFIX" > /dev/null 2> "$WORK/stageG.err" || refusal_status=$?
[[ $refusal_status -eq 1 ]] \
    || fail "G: incompatible options returned status $refusal_status, expected 1" "$WORK/stageG.err"
grep -q -- "--huge-pages cannot be combined with the compact checkpoint" "$WORK/stageG.err" \
    || fail "G: no refusal message" "$WORK/stageG.err"
require_empty_registry

# Page availability can change after either cleanup: classify each attempt on
# its own outcome, never use an earlier success as proof that resources remain.
control_status=0
"$BIN" shm --huge-pages "$PREFIX" > /dev/null 2> "$WORK/stageGcontrol.err" || control_status=$?
if [[ $control_status -eq 0 ]]; then
    require_default_layout
    shm_test_destroy "G: clear the control stage"
    require_empty_registry
elif huge_unavailable "$control_status" "$WORK/stageGcontrol.err"; then
    require_empty_registry
    echo "SKIP: shm_compact_occ_hugepage_test: unset-variable stage unavailable (status $control_status); the --compact-occ refusal passed"
    cat "$WORK/stageGcontrol.err"
else
    fail "G: huge-page control failed (status $control_status)" "$WORK/stageGcontrol.err"
fi

env_status=0
BWA3_COMPACT_OCC=1 "$BIN" shm --huge-pages "$PREFIX" \
    > /dev/null 2> "$WORK/stageG2.err" || env_status=$?
grep -q "compact" "$WORK/stageG2.err" \
    && fail "G: BWA3_COMPACT_OCC=1 changed shm --huge-pages" "$WORK/stageG2.err"
if [[ $env_status -eq 0 ]]; then
    require_default_layout
elif huge_unavailable "$env_status" "$WORK/stageG2.err"; then
    require_empty_registry
    echo "SKIP: shm_compact_occ_hugepage_test: env-enabled stage unavailable (status $env_status); the --compact-occ refusal passed"
    cat "$WORK/stageG2.err"
else
    fail "G: env-enabled stage failed (status $env_status)" "$WORK/stageG2.err"
fi
shm_test_destroy "G: final cleanup"

echo "PASS: shm_compact_occ_hugepage_test"
