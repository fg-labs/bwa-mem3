#!/usr/bin/env bash
# Every standalone and doctest binary must relink when a linked archive changes.
# Requires a matching-flags vendored build: make -C test all. MAKE selects make.
# Dry runs use -W; only the executed link touches timestamps, restored on exit.
# System mode is dependency-tested without requiring an installed system zlib-ng.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."
MAKE="${MAKE:-make}"
libs=(libbwa.a ext/zlib-ng/build/libz-ng.a)
# Parse-time $(info) notes (e.g. the ARM_TUNE fallback) share stdout with
# print-%, so keep only the printed values: the last line of each capture.
if ! bin_list=$("$MAKE" -s -C test print-EXE print-DOCTEST_BINS); then
    echo 'FAIL: cannot enumerate all EXE and DOCTEST_BINS from test/Makefile' >&2
    exit 1
fi
read -r -a bins <<< "$(tail -n 2 <<< "$bin_list" | tr '\n' ' ')"
if ((${#bins[@]} == 0)); then
    echo 'FAIL: test/Makefile has no test binaries'
    exit 1
fi
for file in "${libs[@]}"; do
    if [[ ! -f $file ]]; then
        echo "FAIL: $file is not built - run this against a vendored build"
        exit 1
    fi
done
for bin in "${bins[@]}"; do
    if [[ ! -x test/$bin ]]; then
        echo "FAIL: test/$bin is not built or executable - run make -C test all with matching build flags first"
        exit 1
    fi
done

root=$PWD
# Fixture timestamps are derived from the newest input (or now, if later) so
# their ordering against real inputs holds on any date. Probe stat/date once:
# GNU and BSD spell mtime-in-seconds and epoch formatting differently.
if stat -c %Y . > /dev/null 2>&1; then
    stat_mtime=(stat -c %Y) # GNU coreutils
else
    stat_mtime=(stat -f %m) # BSD / macOS
fi
if date -d @0 +%s > /dev/null 2>&1; then
    touch_stamp() { date -d "@$1" +%Y%m%d%H%M.%S; }
else
    touch_stamp() { date -r "$1" +%Y%m%d%H%M.%S; }
fi
newest=$({
    date +%s
    find src test ext/htslib "${libs[@]}" -type f -exec "${stat_mtime[@]}" {} +
} | sort -n | tail -n 1)
# Ordered stamps, one minute apart, all newer than every current input:
# consumers < in-parent archive rebuild < synthetic link output.
consumer_stamp=$(touch_stamp $((newest + 60)))
archive_stamp=$(touch_stamp $((newest + 120)))
link_stamp=$(touch_stamp $((newest + 180)))
export LINK_STAMP=$link_stamp
refdir="$(mktemp -d)"
for i in "${!libs[@]}"; do touch -r "${libs[$i]}" "$refdir/lib$i"; done
for bin in "${bins[@]}"; do touch -r "test/$bin" "$refdir/$bin"; done
restore_mtimes() {
    local i bin
    for i in "${!libs[@]}"; do touch -r "$refdir/lib$i" "${libs[$i]}"; done
    for bin in "${bins[@]}"; do touch -r "$refdir/$bin" "test/$bin"; done
    rm -rf "$refdir"
}
trap restore_mtimes EXIT

# A mode switch must not dirty the test flag stamp or rebuild every object:
# this test isolates archive dependencies, not the mode-switch rebuild gate.
flags_sig=$("$MAKE" -s -C test print-TEST_FLAGS_SIG | tail -n 1)
# Make the already-built consumers newer than all current inputs. -W below
# simulates archive changes without touching archives or waiting for the clock.
for bin in "${bins[@]}"; do touch "test/$bin"; done
would_link() { # <mode> <binary> [what-if archive]
    local extra=()
    [[ -n ${3:-} ]] && extra=(-W "../$3")
    # The existing ARM rowalt phony gate looks remade under -n; isolate it.
    if ! "$MAKE" -n -C test -o libsais-objs -o ../src/bandedSWA.rowalt.o "${extra[@]}" "ZLIBNG_USE_VENDORED=$1" "TEST_FLAGS_SIG=$flags_sig" "$2" > "$refdir/dry-run" 2>&1; then
        cat "$refdir/dry-run" >&2
        echo 'FAIL: standalone relink dry run failed' >&2
        exit 1
    fi
    grep -qE -- "(^| )-o $2( |\$)" "$refdir/dry-run"
}
for mode in 1 0; do
    for bin in "${bins[@]}"; do
        if would_link "$mode" "$bin"; then
            echo "FAIL: $bin already stale before archive simulation (vendored=$mode)"
            exit 1
        fi
    done
done
for lib in "${libs[@]}"; do
    for mode in 1 0; do
        want=1
        [[ $lib != libbwa.a && $mode == 0 ]] && want=0
        for bin in "${bins[@]}"; do
            got=0
            would_link "$mode" "$bin" "$lib" && got=1
            if ((got != want)); then
                echo "FAIL: $bin relink=$got, expected $want after only $lib changed (vendored=$mode)"
                exit 1
            fi
        done
        echo "checked: ${#bins[@]} binaries after only $lib changed (vendored=$mode, expected relink=$want)"
    done
done

# An archive rebuilt by the ordered parent make must be re-statted in THIS
# invocation, not one invocation later. Use the actual test Makefile, read-only
# existing inputs, a tiny parent builder, and a linker that records its outputs.
fixture="$refdir/in-parent"
mkdir -p "$fixture/test" "$fixture/ext/zlib-ng/build" "$fixture/ext/libsais/src"
cp test/Makefile "$fixture/test/Makefile"
ln -s "$root/src" "$fixture/src"
ln -s "$root/libbwa.a" "$fixture/libbwa.a"
ln -s "$root/ext/htslib" "$fixture/ext/htslib"
for path in "$root"/test/* "$root"/test/.[!.]*; do
    name=${path##*/}
    [[ $name == Makefile || $name == .test-build-flags || ! -e $path ]] && continue
    skip=0
    for bin in "${bins[@]}"; do [[ $name == "$bin" ]] && skip=1; done
    ((skip)) || ln -s "$path" "$fixture/test/$name"
done
cat > "$fixture/Makefile" << 'EOF'
.PHONY: parent ext/libsais/src/libsais.o ext/libsais/src/libsais64.o ext/zlib-ng/build/libz-ng.a src/bandedSWA.rowalt.o
parent:
	@echo parent >> parent-calls
ext/libsais/src/libsais.o ext/libsais/src/libsais64.o src/bandedSWA.rowalt.o: parent ;
ext/zlib-ng/build/libz-ng.a: parent
	@echo archive >> archive-calls
	@if test "$(REBUILD_ARCHIVE)" = 1; then touch -t $(ARCHIVE_STAMP) $@; fi
EOF
cat > "$fixture/link" << 'EOF'
#!/usr/bin/env bash
set -euo pipefail
while (($#)); do
    if [[ $1 == -o ]]; then
        case $2 in *.o) echo "FAIL: fixture attempted to compile $2" >&2; exit 1 ;; esac
        echo "$2" >> ../linked
        touch -t "$LINK_STAMP" "$2"
        exit 0
    fi
    shift
done
echo 'FAIL: synthetic linker invoked without an output' >&2
exit 1
EOF
chmod +x "$fixture/link"
hold_inputs=()
for object in test/*.o test/framework/*.o test/unit/*.o test/integration/*.o test/framework/*.a; do
    [[ -e $object ]] && hold_inputs+=(-o "${object#test/}")
done
for mode in 1 0; do
    for changed in 0 1; do
        printf '%s\n' "$flags_sig" > "$fixture/test/.test-build-flags"
        touch -t 202001010000 "$fixture/test/.test-build-flags"
        touch -t "$(touch_stamp "$newest")" "$fixture/ext/zlib-ng/build/libz-ng.a"
        for bin in "${bins[@]}"; do touch -t "$consumer_stamp" "$fixture/test/$bin"; done
        : > "$fixture/linked"
        : > "$fixture/parent-calls"
        : > "$fixture/archive-calls"
        rc=0
        "$MAKE" -s -C "$fixture/test" "${hold_inputs[@]}" "${bins[@]}" "CXX=$fixture/link" ARM_TUNE= UNIT_EXTRA_OBJS= \
            "TEST_FLAGS_SIG=$flags_sig" "ZLIBNG_USE_VENDORED=$mode" "REBUILD_ARCHIVE=$changed" "ARCHIVE_STAMP=$archive_stamp" > "$refdir/in-parent.log" 2>&1 || rc=$?
        if ((rc != 0)); then
            cat "$refdir/in-parent.log" >&2
            echo "FAIL: in-parent archive fixture failed (exit $rc)" >&2
            exit "$rc"
        fi
        want=0
        ((mode && changed)) && want=${#bins[@]}
        got=$(wc -l < "$fixture/linked" | tr -d ' ')
        parents=$(wc -l < "$fixture/parent-calls" | tr -d ' ')
        archives=$(wc -l < "$fixture/archive-calls" | tr -d ' ')
        if ((got != want || parents != 1 || archives != mode)); then
            cat "$refdir/in-parent.log" "$fixture/linked" >&2
            echo "FAIL: in-parent vendored=$mode changed=$changed: links=$got (want $want), parents=$parents (want 1), archive builds=$archives (want $mode)" >&2
            exit 1
        fi
        if ((want)); then
            for bin in "${bins[@]}"; do
                if ! grep -qxF "$bin" "$fixture/linked"; then
                    echo "FAIL: in-parent rebuild did not relink $bin" >&2
                    exit 1
                fi
            done
        fi
        echo "ok: in-parent vendored=$mode changed=$changed: $got links, one ordered parent"
    done
done

# Execute a real standalone relink as well, preserving its status and logs on
# failure before the EXIT trap removes the scratch directory.
touch -t "$archive_stamp" ext/zlib-ng/build/libz-ng.a
rc=0
"$MAKE" -C test ZLIBNG_USE_VENDORED=1 "TEST_FLAGS_SIG=$flags_sig" packed_text_test > "$refdir/real-link" 2>&1 || rc=$?
if ((rc != 0)); then
    cat "$refdir/real-link" >&2
    echo "FAIL: standalone archive-only relink failed (exit $rc)" >&2
    exit "$rc"
fi
if ! grep -qE -- '(^| )-o packed_text_test( |$)' "$refdir/real-link"; then
    cat "$refdir/real-link" >&2
    echo 'FAIL: standalone make did not execute the archive-only relink'
    exit 1
fi
rc=0
./test/packed_text_test > "$refdir/packed-text" 2>&1 || rc=$?
if ((rc != 0)); then
    cat "$refdir/packed-text" >&2
    echo "FAIL: packed_text_test failed after the archive-only relink (exit $rc)" >&2
    exit "$rc"
fi
echo 'ok: standalone make executed archive-only relink; packed_text_test ran successfully'
echo "PASS: make_test_relink (${#bins[@]} legacy/doctest binaries; external and in-parent archive updates; unchanged/system controls)"
