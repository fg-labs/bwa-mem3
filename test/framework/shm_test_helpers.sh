#!/usr/bin/env bash
# Sourced test helpers. Callers own registry preflight and EXIT trap selection.
# BIN and WORK identify the test executable and private diagnostic directory.

# Required destroys preserve their status and stderr; an EXIT retry uses a
# different file so it cannot overwrite the original failure diagnostic.
shm_test_destroy() {
    local context=$1 status=0
    "$BIN" shm -d > /dev/null 2> "$WORK/clear.err" || status=$?
    if [[ $status -ne 0 ]]; then
        # Even an unreadable diagnostic or failed write must retain the status.
        {
            echo "FAIL: $context: shm -d failed (status $status)"
            cat "$WORK/clear.err"
        } >&2 || true
        exit "$status"
    fi
}

# Best-effort EXIT cleanup must never replace the status that triggered it.
shm_test_exit_cleanup() {
    local status=$1 retry_status=0
    "$BIN" shm -d > /dev/null 2> "$WORK/exit-clear.err" || retry_status=$?
    if [[ $retry_status -ne 0 ]]; then
        {
            echo "EXIT cleanup: shm -d retry failed (status $retry_status); original status $status"
            cat "$WORK/exit-clear.err"
        } >&2 || true
    fi
    shm_test_temp_cleanup "$status"
}

# Temp-only EXIT cleanup is also safe after refusing a foreign registry row.
shm_test_temp_cleanup() {
    local status=$1
    if ! rm -rf "$WORK"; then
        echo "EXIT cleanup: cannot remove temporary directory '$WORK'; original status $status" >&2 || true
    fi
    exit "$status"
}

# Only the stager's documented host/resource refusals permit a skip. The CLI
# returns 1 for these; signals and other errors must not count as unavailable.
huge_unavailable() {
    [[ $1 -eq 1 ]] || return 1
    awk '
        /^\[E::bwa_shm_stage\] --huge-pages needs Linux$/ ||
        /^\[E::bwa_shm_stage\] --huge-pages is refused on AMD x86 hosts: / ||
        /^\[E::bwa_shm_stage\] --huge-pages cannot read kernel.shm_rmid_forced \(\/proc\/sys\/kernel\/shm_rmid_forced\): / ||
        /^\[E::bwa_shm_stage\] --huge-pages cannot work with kernel.shm_rmid_forced=1: / ||
        /^\[E::bwa_shm_stage\] --huge-pages needs \/proc\/self\/ns\/(ipc|pid), to record / ||
        (/^\[E::bwa_shm_stage\] cannot create a System V segment of [0-9]+ x 1 GB pages for / &&
         /\((the 1 GB pool has too few unreserved pages .*|it needs CAP_IPC_LOCK or membership of the group .*|this kernel or CPU has no 1 GB huge page size|the segment exceeds kernel.shmmax, or this host has no 1 GB page size|kernel.shmall or kernel.shmmni is exhausted)\)\. Nothing was staged; rerun without --huge-pages to stage in \/dev\/shm$/) ||
        (/^\[E::bwa_shm_stage\] cannot attach a System V segment of [0-9]+ x 1 GB pages for / &&
         /\(the address space or mapping limit of this process\)\. Nothing was staged; rerun without --huge-pages to stage in \/dev\/shm$/) ||
        (/^\[E::bwa_shm_stage\] cannot fault in the [0-9]+ x 1 GB pages of / &&
         /\(a hugetlb cgroup limit, or a NUMA memory policy that excludes the nodes holding the pool.s free pages\)\. Nothing was staged$/) { known++; next }
        /^\[E::main_shm\] failed to stage .* in shared memory$/ { next }
        /^\[(M|W)::/ { next }
        NF { unknown = 1 }
        END { exit !(known && !unknown) }
    ' "$2"
}
