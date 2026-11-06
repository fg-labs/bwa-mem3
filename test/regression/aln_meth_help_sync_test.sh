#!/usr/bin/env bash
# Controls for the live-help gate: stale snippets cannot hide a broken executable.
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
cat > "$tmp/bwa" << 'STUB'
#!/usr/bin/env bash
if [ "$2" = --help ]; then
    case ${CONTROL:-good} in
        failure) exit 7 ;;
        empty) exit 0 ;;
    esac
    printf '   --meth[=CHEM] shared\n   --meth-context all|cpg\n   --meth-scoring collapsed|genomic|neutral\n   --meth-tags SPEC\n'
    if [ "${CONTROL:-good}" != incomplete ]; then printf '   --set-as-failed f|r\n                 QC-fail\n'; fi
    if [ "${CONTROL:-good}" = drift ] && [ "$1" = aln ]; then printf '                 drift\n'; fi
    printf '   Under --meth, command-specific\n'
    exit 0
fi
case ${CONTROL:-good} in
    refusal-success) printf 'ERROR: refused\n' >&2; exit 0 ;;
    refusal-empty) exit 1 ;;
esac
cpg=0 taps=0 collapsed=0 invalid=0
for arg in "$@"; do
    case $arg in
        *bogus*|*nope*|*XQ*|--set-as-failed=x) invalid=1 ;;
        --meth-context=cpg) cpg=1 ;;
        --meth=taps) taps=1 ;;
        --meth-scoring=collapsed) collapsed=1 ;;
    esac
done
if [ "$invalid" = 1 ] || { [ "$cpg" = 1 ] && { [ "$taps" = 0 ] || [ "$collapsed" = 1 ]; }; }; then
    if [ "${CONTROL:-good}" = refusal-drift ] && [ "$1" = mem ]; then
        printf 'ERROR: different\n' >&2
    else printf 'ERROR: refused\n' >&2; fi
    exit 1
fi
if [ "${CONTROL:-good}" = accepted-rejected ]; then printf 'ERROR: invalid option\n' >&2; exit 1; fi
if [ "$1" = aln ]; then diagnostic='original reference index file'; else diagnostic='seed index'; fi
printf 'ERROR: %s missing; bwa-mem3 index --meth\n' "$diagnostic" >&2
exit 1
STUB
chmod +x "$tmp/bwa"
CONTROL=good BWA_MEM3="$tmp/bwa" bash "$root/test/regression/aln_meth_help_sync.sh" > "$tmp/good.log" 2>&1 || {
    printf 'FAIL: gate rejected good control\n' >&2
    exit 1
}
count=0
for control in failure empty incomplete drift refusal-success refusal-empty refusal-drift accepted-rejected; do
    rc=0
    CONTROL=$control BWA_MEM3="$tmp/bwa" bash "$root/test/regression/aln_meth_help_sync.sh" > "$tmp/$control.log" 2>&1 || rc=$?
    if [ "$rc" = 0 ]; then
        printf 'FAIL: gate accepted %s control\n' "$control" >&2
        exit 1
    fi
    count=$((count + 1))
done
printf 'PASS: live-help gate accepts good control and rejects %s broken controls\n' "$count"
