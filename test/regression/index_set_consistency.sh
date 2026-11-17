#!/usr/bin/env bash
# Structural index sanity checks: lengths, sizes and SA-rate tags, not provenance.
# Equal-length content substitutions remain outside this contract.
# Env: BWA_MEM3 (required executable). All fixtures are generated in private scratch.
set -euo pipefail
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
[[ -n ${BWA_MEM3:-} ]] || fail 'BWA_MEM3 must be set'
case "$BWA_MEM3" in
    /*) BIN=$BWA_MEM3 ;;
    *) BIN="$PWD/$BWA_MEM3" ;;
esac
[[ -x $BIN ]] || fail 'BWA_MEM3 is not executable'
command -v python3 > /dev/null 2>&1 || fail 'python3 is required to generate index fixtures'
# shellcheck source=test/framework/shm_registry_guard.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/shm_registry_guard.sh"
WORK=$(mktemp -d)
BASE="x${WORK##*/}.fa"
cleanup() {
    if shm_registry_holds_only "$BIN" "$BASE" "$BASE.meth"; then
        "$BIN" shm -d > /dev/null 2>&1 || true
    fi
    rm -rf "$WORK"
}
trap cleanup EXIT
cd "$WORK"
python3 - << 'PY'
import random
from pathlib import Path
rng = random.Random(31)
a = ''.join(rng.choice('ACGT') for _ in range(30001))
b = ''.join(rng.choice('ACGT') for _ in range(17003))
Path('a.fa').write_text('>chrA\n' + a + '\n')
Path('b.fa').write_text('>chrB\n' + b + '\n')
for remainder in range(4):
    Path('r%d.fa' % remainder).write_text('>chrR\n' + a[:12000 + remainder] + '\n')
Path('empty-name.fa').write_text('>\n' + a[:1200] + '\n')
Path('long-name.fa').write_text('>' + 'x' * 9000 + '\n' + a[:1200] + '\n')
with open('r.fq', 'w') as out:
    for i in range(50):
        p = rng.randint(0, len(a) - 150)
        out.write('@r%d\n%s\n+\n%s\n' % (i, a[p:p+150], 'I' * 150))
PY
index() { "$BIN" index -t 1 "$@" > build.out 2> build.err || {
    cat build.err >&2
    fail "index $*"
}; }
records() { grep -v '^@PG' "$1" > "$2"; }
index a.fa
index b.fa
"$BIN" mem -t 1 a.fa r.fq > good.raw 2> good.err || fail 'valid disk load'
records good.raw good.sam
[[ $(grep -cv '^@' good.sam) -eq 50 ]] || fail 'valid load wrote fewer than 50 records'
mkset() {
    mkdir -p "$1"
    cp "$2.bwt.2bit.64" "$1/$BASE.bwt.2bit.64"
    local suffix
    for suffix in pac ann amb; do cp "$3.$suffix" "$1/$BASE.$suffix"; done
}
refused() {
    local label=$1 expected=$2
    shift 2
    if "$BIN" mem -t 1 "$@" > "$label.sam" 2> "$label.err"; then fail "$label accepted"; fi
    grep -qF "$expected" "$label.err" || {
        tail -8 "$label.err" >&2
        fail "$label diagnostic"
    }
    [[ $(grep -cv '^@' "$label.sam" || true) -eq 0 ]] || fail "$label wrote records"
    echo "PASS: $label refused"
}
shm_refused() {
    local label=$1 expected=$2
    if "$BIN" shm "$label/$BASE" > "shm-$label.out" 2> "shm-$label.err"; then fail "shm $label accepted"; fi
    grep -qF "$expected" "shm-$label.err" || {
        tail -8 "shm-$label.err" >&2
        fail "shm $label diagnostic"
    }
    local listing
    listing=$("$BIN" shm -l) || fail 'shm listing failed'
    if grep -qF "$BASE" <<< "$listing"; then fail "shm $label staged a rejected set"; fi
    echo "PASS: shm $label refused"
}
for spec in 'long:a.fa:b.fa' 'short:b.fa:a.fa'; do
    IFS=: read -r label bwt metadata <<< "$spec"
    mkset "$label" "$bwt" "$metadata"
    refused "$label" 'do not belong together' "$label/$BASE" r.fq
    refused "compact-$label" 'do not belong together' --compact-occ "$label/$BASE" r.fq
    shm_refused "$label" 'do not belong together'
    before=$(cksum "$label/$BASE.bwt.2bit.64")
    if "$BIN" re-sa -u 2 "$label/$BASE" > resa.out 2> resa.err; then fail "re-sa $label accepted"; fi
    grep -qF 'do not belong together' resa.err || fail "re-sa $label diagnostic"
    [[ $(cksum "$label/$BASE.bwt.2bit.64") == "$before" ]] || fail 're-sa modified refused input'
done
for kind in short long tag; do
    mkset "bwt-$kind" a.fa a.fa
    python3 - "$kind" "bwt-$kind/$BASE.bwt.2bit.64" << 'PY'
import os
import struct
import sys
mode, path = sys.argv[1:]
with open(path, 'r+b') as out:
    if mode == 'short':
        out.truncate(os.path.getsize(path) - 1)
    elif mode == 'long':
        out.seek(0, 2)
        out.write(b'X')
    else:
        out.seek(-8, 2)
        out.write(struct.pack('=q', 99))
PY
    refused "bwt-$kind" 'invalid size, header or SA-rate tag' "bwt-$kind/$BASE" r.fq
    shm_refused "bwt-$kind" 'invalid size, header or SA-rate tag'
    before=$(cksum "bwt-$kind/$BASE.bwt.2bit.64")
    for args in inspect densify; do
        rc=0
        if [[ $args == inspect ]]; then
            "$BIN" re-sa "bwt-$kind/$BASE" > resa.out 2> resa.err || rc=$?
        else
            "$BIN" re-sa -u 2 "bwt-$kind/$BASE" > resa.out 2> resa.err || rc=$?
        fi
        [[ $rc -ne 0 ]] || fail "re-sa accepted bad $kind $args"
        grep -qF 'invalid size, header or SA-rate tag' resa.err || fail 're-sa size diagnostic'
        [[ $(cksum "bwt-$kind/$BASE.bwt.2bit.64") == "$before" ]] || fail 're-sa changed invalid file'
    done
done
# Check every rounding residue and a same-size .pac with a changed trailing byte.
for remainder in 0 1 2 3; do
    index "r$remainder.fa"
    "$BIN" mem -t 1 "r$remainder.fa" r.fq > /dev/null 2> /dev/null || fail 'valid remainder refused'
    mkset "pac-$remainder" "r$remainder.fa" "r$remainder.fa"
    python3 - "pac-$remainder/$BASE.pac" << 'PY'
import sys
with open(sys.argv[1], 'r+b') as out:
    out.seek(-1, 2)
    value = out.read(1)[0]
    out.seek(-1, 2)
    out.write(bytes([(value + 1) % 4]))
PY
    refused "pac-$remainder" 'invalid length remainder' "pac-$remainder/$BASE" r.fq
    shm_refused "pac-$remainder" 'invalid length remainder'
done
for kind in short long; do
    mkset "pac-$kind" a.fa a.fa
    python3 - "$kind" "pac-$kind/$BASE.pac" << 'PY'
import os
import sys
mode, path = sys.argv[1:]
with open(path, 'r+b') as out:
    out.truncate(os.path.getsize(path) + (1 if mode == 'long' else -1))
PY
    refused "pac-$kind" "$BASE.pac is" "pac-$kind/$BASE" r.fq
    shm_refused "pac-$kind" "$BASE.pac is"
done

# Invalid packed data must be rejected before compact record conversion.
mkset pac-first a.fa a.fa
printf 'X' >> "pac-first/$BASE.pac"
python3 - "pac-first/$BASE.bwt.2bit.64" << 'PYTHON'
import struct
import sys
with open(sys.argv[1], 'r+b') as out:
    out.seek(48)
    out.write(struct.pack('=q', -1))
PYTHON
refused pac-first "$BASE.pac is" --compact-occ "pac-first/$BASE" r.fq

# Strip only the rate-3 tag: this reproduces the stock bwa-mem2 disk layout.
mkset legacy a.fa a.fa
python3 - "legacy/$BASE.bwt.2bit.64" << 'PY'
import os
import sys
os.truncate(sys.argv[1], os.path.getsize(sys.argv[1]) - 8)
PY
"$BIN" mem -t 1 "legacy/$BASE" r.fq > legacy.raw 2> legacy.err || fail 'legacy untagged disk load'
records legacy.raw legacy.sam
cmp good.sam legacy.sam || fail 'legacy output differs'
"$BIN" mem --compact-occ -t 1 "legacy/$BASE" r.fq > legacy.raw 2> legacy.err || fail 'legacy compact load'
records legacy.raw legacy.sam
cmp good.sam legacy.sam || fail 'legacy compact output differs'
"$BIN" re-sa "legacy/$BASE" > /dev/null 2> legacy-resa.err || fail 'legacy re-sa inspection'
"$BIN" re-sa -u 2 "legacy/$BASE" > /dev/null 2> legacy-resa.err || fail 'legacy re-sa densify'
"$BIN" mem -t 1 "legacy/$BASE" r.fq > legacy.raw 2> legacy.err || fail 'resampled load'
records legacy.raw legacy.sam
cmp good.sam legacy.sam || fail 'resampled output differs'
echo 'PASS: legacy bwa-mem2 layout and re-sa output accepted'

# All successful global-shm operations are isolated behind the registry guard.
shm_registry_guard "$BIN" index_set_consistency "$BASE" "$BASE.meth"
"$BIN" shm "legacy/$BASE" > /dev/null 2> legacy-shm.err || fail 'legacy/resampled shm stage'
"$BIN" mem -t 1 "legacy/$BASE" r.fq > legacy.raw 2> legacy.err || fail 'resampled shm attach'
records legacy.raw legacy.sam
cmp good.sam legacy.sam || fail 'attached output differs'
"$BIN" shm -d > /dev/null 2>&1 || fail 'drop own stage'
# Stage the genuinely untagged file too, not just the resampled tagged one.
mkset legacy-stage a.fa a.fa
python3 - "legacy-stage/$BASE.bwt.2bit.64" << 'PY'
import os
import sys
os.truncate(sys.argv[1], os.path.getsize(sys.argv[1]) - 8)
PY
"$BIN" shm "legacy-stage/$BASE" > /dev/null 2> legacy-shm.err || fail 'untagged shm stage'
"$BIN" mem -t 1 "legacy-stage/$BASE" r.fq > legacy.raw 2> legacy.err || fail 'untagged shm attach'
records legacy.raw legacy.sam
cmp good.sam legacy.sam || fail 'legacy shm output differs'
"$BIN" shm -d > /dev/null 2>&1 || fail 'drop own untagged stage'

mkdir meth
cp a.fa "meth/$BASE"
index --meth "meth/$BASE"
"$BIN" shm --meth "meth/$BASE" > /dev/null 2> meth-shm.err || fail 'valid meth stage'
# Replacing only the original leaves seed metadata from a different reference.
index --force -p "meth/$BASE" b.fa
refused meth-staged 'inconsistent reference length' --meth "meth/$BASE" r.fq
grep -qF 'shm -d' meth-staged.err || fail 'staged mismatch omitted shm -d remedy'
grep -qF 'index --meth --force' meth-staged.err || fail 'meth mismatch omitted --meth remedy'
"$BIN" shm -d > /dev/null 2>&1 || fail 'drop own meth stage'
refused meth-disk 'inconsistent reference length' --meth "meth/$BASE" r.fq
if "$BIN" shm --meth "meth/$BASE" > /dev/null 2> meth-bad-shm.err; then fail 'staged mismatched meth metadata'; fi
grep -qF 'inconsistent reference length' meth-bad-shm.err || fail 'meth stage mismatch diagnostic'
for name in empty-name long-name; do
    index --meth "$name.fa"
    "$BIN" mem --meth -t 1 "$name.fa" r.fq > "$name.sam" 2> "$name.err" || {
        tail -8 "$name.err" >&2
        fail "valid $name meth load"
    }
done
echo 'PASS: meth metadata mismatch refused, valid upstream name round trips accepted'
echo 'PASS: index_set_consistency'
