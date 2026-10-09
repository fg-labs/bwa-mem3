#!/usr/bin/env bash
# Source-only helper. Resolve an external timer with the platform's RSS output.
# Sets TIME_CMD; the optional fallback argument is for controlled regression tests.
bwamem3_resolve_time() {
    local platform="$1" fallback="${2-/usr/bin/time}" flag candidate report
    local path_time
    path_time="$(type -P time || true)"
    case "$platform" in
        Darwin) flag=-l ;;
        *) flag=-v ;;
    esac
    TIME_CMD=()
    for candidate in "$path_time" "$fallback"; do
        [[ -n "$candidate" && -f "$candidate" && -x "$candidate" ]] || continue
        if ! report="$(LC_ALL=C "$candidate" "$flag" /bin/sh -c ':' 2>&1)"; then
            continue
        fi
        if [[ "$platform" == Darwin ]]; then
            [[ "$report" =~ [0-9]+[[:space:]]+maximum[[:space:]]resident[[:space:]]set[[:space:]]size ]] || continue
            [[ "$report" =~ [0-9.]+[[:space:]]+real ]] || continue
        else
            [[ "$report" =~ Maximum[[:space:]]resident[[:space:]]set[[:space:]]size[[:space:]]\(kbytes\):[[:space:]]*[0-9]+ ]] || continue
            [[ "$report" == *"Elapsed (wall clock) time"* ]] || continue
        fi
        # Caller consumes the resolved command array.
        # shellcheck disable=SC2034
        TIME_CMD=("$candidate" "$flag")
        return 0
    done
    printf 'FAIL: no compatible external time command; install GNU time on Linux or BSD time on macOS and add it to PATH.\n' >&2
    return 1
}
