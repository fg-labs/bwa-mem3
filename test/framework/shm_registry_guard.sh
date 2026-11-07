#!/usr/bin/env bash
# test/framework/shm_registry_guard.sh
#
# Source-only helper for tests that run the global `bwa-mem3 shm -d`, which
# drops every index staged on the host, not just the test's own. Before its
# first `shm -d`, such a test (which has already done `cd "$(dirname "$0")/.."`,
# so it runs from the repository root) calls
#
#     # shellcheck source=test/framework/shm_registry_guard.sh
#     source test/framework/shm_registry_guard.sh
#     shm_registry_guard "$BIN" <test-name> <basename>...
#
# which returns when the registry is empty or lists only the given basenames
# (what the test itself stages, possibly left by an earlier aborted run), and
# otherwise prints SKIP and exits 0. It reads the registry with `shm -l` only,
# so the check itself changes nothing. A failing `shm -l` counts as a foreign
# registry. The C++ tests use test/framework/shm_registry_guard.h.
#
# Portability: bash 3.2 (macOS /bin/bash); no mapfile or associative arrays.

# shm_registry_holds_only <bin> <basename>...: status 0 when `<bin> shm -l`
# succeeds and lists no index other than the given basenames.
shm_registry_holds_only() {
    local bin=$1 listing name allowed ok
    shift
    listing=$("$bin" shm -l 2> /dev/null) || return 1
    while IFS=$'\t' read -r name _; do
        [[ -z "$name" ]] && continue
        ok=0
        for allowed in "$@"; do
            if [[ "$name" == "$allowed" ]]; then
                ok=1
                break
            fi
        done
        [[ "$ok" -eq 1 ]] || return 1
    done <<< "$listing"
    return 0
}

# shm_registry_guard <bin> <test-name> <basename>...: SKIP (exit 0) unless the
# registry holds only the given basenames.
shm_registry_guard() {
    local bin=$1 test_name=$2
    shift 2
    if ! shm_registry_holds_only "$bin" "$@"; then
        echo "SKIP: $test_name: the shared-memory registry lists an index other than" \
            "'$*', or \`$bin shm -l\` failed. This test runs \`shm -d\`, which would drop" \
            "every staged index on the host. Inspect with \`bwa-mem3 shm -l\`; clear with" \
            "\`bwa-mem3 shm -d\` or run on another host." >&2
        echo "SKIP: $test_name (shared-memory registry holds other indices)"
        exit 0
    fi
}
