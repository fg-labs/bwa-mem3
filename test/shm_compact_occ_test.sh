#!/usr/bin/env bash
# End-to-end regression for the compact checkpoint layout in shared memory
# (`shm --compact-occ`), on the phix index. A staged index is attached in
# whichever layout it holds; --compact-occ and BWA3_COMPACT_OCC choose the
# layout of a disk load only.
#   A. a segment staged with --compact-occ is listed as `compact`, a
#      `mem --compact-occ` run attaches it quietly, and its SAM is
#      byte-identical to a disk-loaded default run (also under
#      --compat=bwa-mem2);
#   B. a run that requested no layout attaches it too, with an [M::] note;
#      BWA3_COMPACT_OCC=1 attaches it quietly, and BWA3_COMPACT_OCC=0 (the
#      default layout, requested) attaches it with a [W::] line;
#   C. staging the default layout over it fails (one layout per index);
#   D. the converse: a default segment is attached quietly by a default run,
#      and with a [W::] line by `mem --compact-occ` and BWA3_COMPACT_OCC=1;
#      `shm --compact-occ` over it fails;
#   E. `shm` does not read BWA3_COMPACT_OCC: with it set, `shm` stages the
#      default layout, and `shm -u` is not refused;
#   F. `shm --compact-occ` with -u (any value), -l or -d is refused;
#   G. huge-page checks live in shm_compact_occ_hugepage_test.sh so a host
#      refusal there does not mark the normal-layout checks here skipped;
#   H. a seed-only `shm --meth --compact-occ` stage, attached by
#      `mem --meth`, gives the disk-loaded run's SAM.
# @PG is excluded from the comparisons (its CL: records argv).

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR/.."

BIN=${BIN:-./bwa-mem3}
PREFIX=test/fixtures/phix.fa
READS=test/fixtures/reads.fa

if [[ ! -x "$BIN" ]]; then
    echo "FAIL: $BIN not built. Run 'make -j4' first." >&2
    exit 2
fi

# `shm -d` drops every index staged on the host. Refuse (SKIP) unless the
# registry is empty or holds only this test's own indices.
# shellcheck source=test/framework/shm_registry_guard.sh
source "$SCRIPT_DIR/framework/shm_registry_guard.sh"
# shellcheck source=test/framework/shm_test_helpers.sh
source "$SCRIPT_DIR/framework/shm_test_helpers.sh"
shm_registry_guard "$BIN" shm_compact_occ_test phix.fa ref.fa.meth

need_index=0
for ext in .amb .ann .bwt.2bit.64 .pac; do
    if [[ ! -s "${PREFIX}${ext}" ]]; then
        need_index=1
        break
    fi
done
if [[ "$need_index" -eq 1 ]]; then
    echo "[setup] Building phix index..."
    # --force: a partial fixture index (some files missing) is rebuilt in full;
    # without it `index` refuses to replace the files that remain.
    "$BIN" index --force "$PREFIX" > /dev/null || {
        echo "FAIL: index failed for '$PREFIX'" >&2
        exit 1
    }
fi

WORK=$(mktemp -d -t shm_compact.XXXXXX)
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
body() { grep -v '^@PG' "$1"; }
# run <tag> <mem args...>: SAM to $WORK/<tag>.sam, stderr to $WORK/<tag>.err
run() {
    local tag=$1
    shift
    "$BIN" mem "$@" "$PREFIX" "$READS" > "$WORK/$tag.sam" 2> "$WORK/$tag.err" \
        || fail "$tag: mem exited non-zero" "$WORK/$tag.err"
}
# run_env <tag> <BWA3_COMPACT_OCC value> <mem args...>
run_env() {
    local tag=$1 val=$2
    shift 2
    BWA3_COMPACT_OCC=$val "$BIN" mem "$@" "$PREFIX" "$READS" > "$WORK/$tag.sam" 2> "$WORK/$tag.err" \
        || fail "$tag: mem exited non-zero" "$WORK/$tag.err"
}
same() {
    local tag
    for tag in "$1" "$2"; do
        awk 'NF && !/^@/ { found = 1 } END { exit !found }' "$WORK/$tag.sam" \
            || fail "$tag: SAM has no alignment records" "$WORK/$tag.err"
    done
    cmp -s <(body "$WORK/$1.sam") <(body "$WORK/$2.sam") \
        || fail "$1 vs $2: SAM records differ" "$WORK/$1.err" "$WORK/$2.err"
}
attached() { grep -q "attached from shm" "$WORK/$1.err"; }
# attached_as <tag> <compact|default>: attached, in that layout
attached_as() {
    attached "$1" || fail "$1: did not attach the staged segment" "$WORK/$1.err"
    if [[ $2 == compact ]]; then
        grep -q "attached from shm.*(compact checkpoint layout)" "$WORK/$1.err" \
            || fail "$1: did not attach in the compact layout" "$WORK/$1.err"
    else
        grep -q "attached from shm.*(compact checkpoint layout)" "$WORK/$1.err" \
            && fail "$1: attached in the compact layout" "$WORK/$1.err"
    fi
    return 0
}
warned() { grep -q '^\[W::load_index\] the .* checkpoint layout was requested' "$WORK/$1.err"; }
noted() { grep -q '^\[M::load_index\] the index staged in shared memory holds the compact' "$WORK/$1.err"; }
quiet() {
    warned "$1" && fail "$1: unexpected layout warning" "$WORK/$1.err"
    noted "$1" && fail "$1: unexpected layout note" "$WORK/$1.err"
    return 0
}
list_registry() {
    "$BIN" shm -l > "$WORK/list.out" 2> "$WORK/list.err" \
        || fail "shm -l failed" "$WORK/list.err" "$WORK/list.out"
}
require_layout() {
    list_registry
    awk -F'\t' -v n="$1" -v layout="$2" '
        $1 == n {
            count++
            if ($2 ~ /^[0-9]+$/ && $2 > 0 &&
                ((layout == "default" && NF == 2) ||
                 (layout == "compact" && NF == 3 && $3 == "compact"))) valid++
        }
        END { exit !(count == 1 && valid == 1) }
    ' "$WORK/list.out" || fail "shm -l: expected exactly one $2 row for $1" "$WORK/list.out"
}
require_empty_registry() {
    list_registry
    [[ ! -s "$WORK/list.out" ]] || fail "shm -l: expected empty registry" "$WORK/list.out"
}

echo "[base] disk-loaded default runs"
run disk
run disk_compat --compat=bwa-mem2
attached disk && fail "disk: attached with nothing staged" "$WORK/disk.err"

echo "[A] shm --compact-occ, then mem --compact-occ attaches it"
"$BIN" shm --compact-occ -t 2 "$PREFIX" > /dev/null 2> "$WORK/stageA.err" \
    || fail "A: shm --compact-occ exited non-zero" "$WORK/stageA.err"
require_layout phix.fa compact
run A --compact-occ
run A_compat --compact-occ --compat=bwa-mem2
for t in A A_compat; do
    attached_as "$t" compact
    quiet "$t"
done
same A disk
same A_compat disk_compat

echo "[B] a compact segment is attached whatever the run requested"
run B
attached_as B compact
noted B || fail "B: no note naming the attached compact layout" "$WORK/B.err"
warned B && fail "B: warned with no layout requested" "$WORK/B.err"
same B disk
run_env B_env1 1
attached_as B_env1 compact
quiet B_env1
same B_env1 disk
run_env B_env0 0
attached_as B_env0 compact
warned B_env0 || fail "B: BWA3_COMPACT_OCC=0 run did not warn" "$WORK/B_env0.err"
same B_env0 disk

echo "[C] staging the default layout over a compact segment fails"
if "$BIN" shm "$PREFIX" > /dev/null 2> "$WORK/stageC.err"; then
    fail "C: shm (default) over a compact segment succeeded" "$WORK/stageC.err"
fi
grep -q "already in shm with the compact" "$WORK/stageC.err" \
    || fail "C: no layout-mismatch error" "$WORK/stageC.err"
"$BIN" shm --compact-occ "$PREFIX" > /dev/null 2> "$WORK/stageC2.err" \
    || fail "C: restaging the same compact layout failed" "$WORK/stageC2.err"

echo "[D] a default segment is attached whatever the run requested"
shm_test_destroy "D"
"$BIN" shm "$PREFIX" > /dev/null 2> "$WORK/stageD.err" \
    || fail "D: shm exited non-zero" "$WORK/stageD.err"
require_layout phix.fa default
run D_def
attached_as D_def default
quiet D_def
same D_def disk
run D --compact-occ
run D_compat --compact-occ --compat=bwa-mem2
run_env D_env1 1
for t in D D_compat D_env1; do
    attached_as "$t" default
    warned "$t" || fail "$t: no warning that the requested compact layout is not used" "$WORK/$t.err"
done
same D disk
same D_compat disk_compat
same D_env1 disk
if "$BIN" shm --compact-occ "$PREFIX" > /dev/null 2> "$WORK/stageD2.err"; then
    fail "D: shm --compact-occ over a default segment succeeded" "$WORK/stageD2.err"
fi
grep -q "already in shm with the default" "$WORK/stageD2.err" \
    || fail "D: no layout-mismatch error" "$WORK/stageD2.err"

echo "[E] shm does not read BWA3_COMPACT_OCC"
shm_test_destroy "E: before env-enabled stage"
BWA3_COMPACT_OCC=1 "$BIN" shm "$PREFIX" > /dev/null 2> "$WORK/stageE.err" \
    || fail "E: BWA3_COMPACT_OCC=1 shm exited non-zero" "$WORK/stageE.err"
require_layout phix.fa default
shm_test_destroy "E: before env-enabled densify"
BWA3_COMPACT_OCC=1 "$BIN" shm -u 2 "$PREFIX" > /dev/null 2> "$WORK/stageE2.err" \
    || fail "E: BWA3_COMPACT_OCC=1 shm -u exited non-zero" "$WORK/stageE2.err"
require_layout phix.fa default
shm_test_destroy "E: after env-enabled densify"

echo "[F] shm --compact-occ with -u, -l or -d is refused"
for u in 2 6; do
    if "$BIN" shm --compact-occ -u "$u" "$PREFIX" > /dev/null 2> "$WORK/stageF$u.err"; then
        fail "F: shm --compact-occ -u $u succeeded" "$WORK/stageF$u.err"
    fi
    grep -q -- "--compact-occ cannot be combined with -u" "$WORK/stageF$u.err" \
        || fail "F: no refusal message for -u $u" "$WORK/stageF$u.err"
    require_empty_registry
done
for flag in -l -d; do
    if "$BIN" shm --compact-occ "$flag" > /dev/null 2> "$WORK/stageF$flag.err"; then
        fail "F: shm --compact-occ $flag succeeded" "$WORK/stageF$flag.err"
    fi
    grep -q -- "--compact-occ cannot be combined with -l or -d" "$WORK/stageF$flag.err" \
        || fail "F: no refusal message for $flag" "$WORK/stageF$flag.err"
done

echo "[H] shm --meth --compact-occ (seed-only), attached by mem --meth"
cp test/fixtures/phix.fa "$WORK/ref.fa"
"$BIN" index --meth "$WORK/ref.fa" > /dev/null || fail "H: index --meth failed for '$WORK/ref.fa'"
"$BIN" mem --meth "$WORK/ref.fa" "$READS" > "$WORK/H_disk.sam" 2> "$WORK/H_disk.err" \
    || fail "H: disk mem --meth exited non-zero" "$WORK/H_disk.err"
"$BIN" shm --meth --compact-occ "$WORK/ref.fa" > /dev/null 2> "$WORK/stageH.err" \
    || fail "H: shm --meth --compact-occ exited non-zero" "$WORK/stageH.err"
require_layout ref.fa.meth compact
"$BIN" mem --meth --compact-occ "$WORK/ref.fa" "$READS" > "$WORK/H.sam" 2> "$WORK/H.err" \
    || fail "H: mem --meth --compact-occ exited non-zero" "$WORK/H.err"
attached_as H compact
quiet H
same H H_disk
shm_test_destroy "H: final cleanup"

echo "PASS: shm_compact_occ_test"
