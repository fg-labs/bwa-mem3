#!/usr/bin/env bash
# Hermetic timer-discovery checks; no build or indexing is needed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TIME_CMD=()
# shellcheck source=scripts/time_command.sh
source "$ROOT/scripts/time_command.sh"
TD="$(mktemp -d)"
trap 'rm -rf "$TD"' EXIT
mkdir -p "$TD/path with spaces" "$TD/fallback" "$TD/empty"
cat > "$TD/path with spaces/time" << 'TIMER'
#!/bin/bash
case "$1" in
    -v)
        echo 'Maximum resident set size (kbytes): 123' >&2
        echo 'Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.01' >&2
        ;;
    -l)
        echo '123 maximum resident set size' >&2
        echo '0.01 real' >&2
        ;;
    *) exit 2 ;;
esac
shift
"$@"
TIMER
chmod +x "$TD/path with spaces/time"
cp "$TD/path with spaces/time" "$TD/fallback/time"
check() {
    local platform="$1" search="$2" fallback="$3" expected="$4" flag="$5"
    PATH="$search" bwamem3_resolve_time "$platform" "$fallback"
    [[ "${TIME_CMD[0]}" == "$expected" && "${TIME_CMD[1]}" == "$flag" ]] || {
        echo 'FAIL: timer selection or platform flag differs' >&2
        exit 1
    }
}
check Linux "$TD/path with spaces" "$TD/fallback/time" "$TD/path with spaces/time" -v
check Darwin "$TD/path with spaces" "$TD/fallback/time" "$TD/path with spaces/time" -l
check Linux "$TD/empty" "$TD/fallback/time" "$TD/fallback/time" -v
# A zero exit without parseable RSS is incompatible too.
printf '#!/bin/bash\nexit 0\n' > "$TD/path with spaces/time"
check Linux "$TD/path with spaces" "$TD/fallback/time" "$TD/fallback/time" -v
printf '#!/bin/bash\nexit 2\n' > "$TD/path with spaces/time"
check Darwin "$TD/path with spaces" "$TD/fallback/time" "$TD/fallback/time" -l
if PATH="$TD/path with spaces" bwamem3_resolve_time Linux "$TD/missing" 2> "$TD/error"; then
    echo 'FAIL: incompatible timer accepted' >&2
    exit 1
fi
if PATH="$TD/empty" bwamem3_resolve_time Linux "$TD/missing" 2> "$TD/error"; then
    echo 'FAIL: missing timer accepted' >&2
    exit 1
fi
grep -q 'install GNU time' "$TD/error"
chmod -x "$TD/path with spaces/time"
check Linux "$TD/path with spaces" "$TD/fallback/time" "$TD/fallback/time" -v
check Linux "$TD/empty" "$TD/fallback/time" "$TD/fallback/time" -v
rc=0
"${TIME_CMD[@]}" /bin/sh -c 'exit 7' > /dev/null 2>&1 || rc=$?
[[ "$rc" == 7 ]] || {
    echo 'FAIL: timer lost command exit status' >&2
    exit 1
}
echo 'PASS: external timer precedence, fallback, compatibility, spaces and platform flags'
