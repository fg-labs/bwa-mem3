#!/usr/bin/env bash
# End-to-end check of the shm huge-page warning on Linux: stage an index whose
# segment is larger than one 2 MB huge page (the 1 MB synthetic fixture) and check
# that `shm` and `mem` print the [W::bwa_shm_stage] / [W::bwa_shm_attach]
# line exactly when it can be on huge pages and this host's /dev/shm mount
# leaves it on small pages (huge=never, or no huge= option; shmem_enabled not
# [deny]; a segment of at least one PMD page), once per process, never at
# -v 1, and that the SAM records do not change. Where smaps cannot show the
# segment's page size, the warning must stay silent. Where the
# host permits a 1 GB stage, `shm --huge-pages` is attached without the warning.
# Documented unavailable hugetlb prerequisites skip that leg with diagnostics.
# SKIP off Linux, where the warning is compiled out.

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR/.."

BIN=${BIN:-./bwa-mem3}
FASTA=test/fixtures/synthetic_1mb.fa

if [[ ! -x "$BIN" ]]; then
    echo "FAIL: $BIN not built. Run 'make -j4' first." >&2
    exit 2
fi
if [[ "$(uname -s)" != Linux ]]; then
    echo "SKIP: shm_hugepage_warn_test: the warning is Linux only"
    exit 0
fi

# shellcheck source=test/framework/shm_registry_guard.sh
source "$SCRIPT_DIR/framework/shm_registry_guard.sh"
# shellcheck source=test/framework/shm_test_helpers.sh
source "$SCRIPT_DIR/framework/shm_test_helpers.sh"
shm_registry_guard "$BIN" shm_hugepage_warn_test hpw.fa

WORK=$(mktemp -d -t shm_hpw.XXXXXX)
trap 'shm_test_exit_cleanup "$?"' EXIT
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
count() { grep -c "^\[W::$1\] the shm segment for 'hpw.fa' is on" "$2" || true; }
same() {
    local sam
    for sam in "$1" "$2"; do
        awk 'NF && !/^@/ { found = 1 } END { exit !found }' "$sam" \
            || fail "$sam: SAM has no alignment records" "$3"
    done
    cmp -s <(grep -v '^@PG' "$1") <(grep -v '^@PG' "$2") \
        || fail "$1 vs $2: SAM records differ" "$3"
}

PREFIX=$WORK/hpw.fa
cp "$FASTA" "$PREFIX"
"$BIN" index "$PREFIX" > /dev/null || fail "index failed for '$PREFIX'"
# Reads: the first 100 bp of each 1 kb of the reference, as FASTA.
awk '!/^>/ { seq = seq $0 } END { for (i = 0; i + 100 <= length(seq) && n < 200; i += 1000) { n++; print ">r" n; print substr(seq, i + 1, 100) } }' \
    "$FASTA" > "$WORK/reads.fa"

"$BIN" mem "$PREFIX" "$WORK/reads.fa" > "$WORK/disk.sam" 2> "$WORK/disk.err" \
    || fail "disk mem exited non-zero" "$WORK/disk.err"

echo "[A] shm: warns after staging exactly when the segment can be on huge pages and is not"
"$BIN" shm "$PREFIX" > /dev/null 2> "$WORK/stage.err" || fail "shm exited non-zero" "$WORK/stage.err"
# The expectation mirrors shm_warn_if_small_pages: no warning for a segment
# smaller than one PMD page (512 MB on an arm64 kernel with 64 KB pages), with
# shmem THP absent or [deny] (never huge), or when /proc/self/smaps is unreadable
# or has no PMD-mapped fields (the page size cannot be told). [force] and a
# huge= mount put the segment on huge pages unless the host is too fragmented
# to, and the warning fires on that fallback, so either count is accepted
# there; huge=never (the default) leaves it on small pages.
opts=$(awk '$2 == "/dev/shm" { print $4 }' /proc/mounts | tail -1)
shmem=$(cat /sys/kernel/mm/transparent_hugepage/shmem_enabled 2> /dev/null || echo none)
pmd=$(cat /sys/kernel/mm/transparent_hugepage/hpage_pmd_size 2> /dev/null || echo 2097152)
size=$("$BIN" shm -l | awk -F'\t' '$1 == "hpw.fa" { print $2 }')
[[ -n "$size" ]] || fail "shm -l does not list hpw.fa"
smaps=none
grep -qsE '^(ShmemPmdMapped|FilePmdMapped):' /proc/self/smaps && smaps=pmd
if [[ $size -lt $pmd || $shmem == none || $shmem == *"[deny]"* || $smaps == none ]]; then
    want=0
elif [[ $shmem == *"[force]"* ]]; then
    want=any
else
    case ",$opts," in
        *,huge=always,* | *,huge=within_size,* | *,huge=advise,*) want=any ;;
        *) want=1 ;;
    esac
fi
echo "[setup] /dev/shm options: ${opts:-?}; shmem_enabled: $shmem; smaps PMD fields: $smaps; PMD page: $pmd; segment: $size; expect warning: $want"
expect() { # expect <count> <what> <file>
    [[ $want == any && $1 -le 1 ]] && return 0
    [[ "$1" == "$want" ]] || fail "$2: expected $want warning(s), got $1" "$3"
}
expect "$(count bwa_shm_stage "$WORK/stage.err")" shm "$WORK/stage.err"

echo "[B] mem: warns once on attach, the same way; never at -v 1; SAM unchanged"
"$BIN" mem "$PREFIX" "$WORK/reads.fa" > "$WORK/att.sam" 2> "$WORK/att.err" \
    || fail "mem exited non-zero" "$WORK/att.err"
grep -q "attached from shm" "$WORK/att.err" || fail "mem did not attach" "$WORK/att.err"
expect "$(count bwa_shm_attach "$WORK/att.err")" mem "$WORK/att.err"
same "$WORK/att.sam" "$WORK/disk.sam" "$WORK/att.err"
"$BIN" mem -v 1 "$PREFIX" "$WORK/reads.fa" > /dev/null 2> "$WORK/quiet.err" \
    || fail "mem -v 1 exited non-zero" "$WORK/quiet.err"
[[ "$(count bwa_shm_attach "$WORK/quiet.err")" == 0 ]] || fail "mem -v 1 warned" "$WORK/quiet.err"
shm_test_destroy "B: clear the small-page stage"

echo "[C] a 1 GB (hugetlb) stage is attached without the warning"
# Free pages alone do not prove that reservations, permissions and host policy
# allow this stage. Try it and classify only the documented refusal diagnostics.
stage_status=0
"$BIN" shm --huge-pages "$PREFIX" > /dev/null 2> "$WORK/hstage.err" || stage_status=$?
if [[ $stage_status -eq 0 ]]; then
    "$BIN" mem "$PREFIX" "$WORK/reads.fa" > "$WORK/huge.sam" 2> "$WORK/huge.err" \
        || fail "mem (1 GB stage) exited non-zero" "$WORK/huge.err"
    grep -q "attached from shm" "$WORK/huge.err" || fail "mem did not attach the 1 GB stage" "$WORK/huge.err"
    [[ "$(count bwa_shm_attach "$WORK/huge.err")" == 0 ]] \
        || fail "mem warned on a 1 GB stage" "$WORK/huge.err"
    same "$WORK/huge.sam" "$WORK/disk.sam" "$WORK/huge.err"
    shm_test_destroy "C: clear the huge-page stage"
elif huge_unavailable "$stage_status" "$WORK/hstage.err"; then
    echo "SKIP: shm_hugepage_warn_test: hugetlb stage unavailable (status $stage_status); legs A and B passed"
    cat "$WORK/hstage.err"
else
    fail "hugetlb stage failed (status $stage_status)" "$WORK/hstage.err"
fi

echo "PASS: shm_hugepage_warn_test"
