#!/usr/bin/env bash
# Compare live shared help and shared-option validation; never trust stale generated snippets.
set -euo pipefail
bin=${BWA_MEM3:-./bwa-mem3}
[[ $bin = /* ]] || bin=$PWD/$bin
[ -x "$bin" ] || {
    printf 'FAIL: binary is not executable: %s\n' "$bin" >&2
    exit 1
}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
for command in mem aln; do
    "$bin" "$command" --help > "$tmp/$command.help" 2>&1 || {
        printf 'FAIL: %s --help failed\n' "$command" >&2
        exit 1
    }
    mawk '
        /^ *--meth\[=CHEM\]/ {on=1; start++}
        on && /^ *--(meth-seed-prune|chimera-qc)/ {exit}
        on && /^ *Under --meth,/ {exit}
        on {
            print
            if (/^ *--meth-context /) context++
            if (/^ *--meth-scoring /) scoring++
            if (/^ *--meth-tags /) tags++
            if (/^ *--set-as-failed f\|r/) failed++
        }
        END {if(start!=1 || context!=1 || scoring!=1 || tags!=1 || failed!=1) exit 1}
    ' "$tmp/$command.help" > "$tmp/$command.block" || {
        printf 'FAIL: missing/incomplete shared block in live %s help\n' "$command" >&2
        exit 1
    }
done
if ! cmp -s "$tmp/mem.block" "$tmp/aln.block"; then
    printf 'FAIL: live shared --meth help differs\n' >&2
    diff -u "$tmp/mem.block" "$tmp/aln.block" >&2
    exit 1
fi
refused=0
for args in '--meth=bogus' '--meth-context=nope' '--meth-scoring=nope' '--meth-tags=XQ' \
    '--meth-context=cpg' '--meth --meth-context=cpg' \
    '--meth=taps --meth-context=cpg --meth-scoring=collapsed' \
    '--meth --set-as-failed=x'; do
    read -r -a options <<< "$args"
    for command in mem aln; do
        rc=0
        "$bin" "$command" "${options[@]}" /nonexistent /nonexistent > "$tmp/$command.out" 2> "$tmp/$command.err" || rc=$?
        [ "$rc" = 1 ] && [ ! -s "$tmp/$command.out" ] || {
            printf 'FAIL: %s %s did not refuse cleanly (exit %s)\n' "$command" "$args" "$rc" >&2
            exit 1
        }
        mawk '/^ERROR/ {print; found=1; exit} END {if(!found) exit 1}' "$tmp/$command.err" > "$tmp/$command.error" || {
            printf 'FAIL: %s %s has no diagnostic\n' "$command" "$args" >&2
            exit 1
        }
    done
    cmp -s "$tmp/mem.error" "$tmp/aln.error" || {
        printf 'FAIL: refusal differs: %s\n' "$args" >&2
        exit 1
    }
    refused=$((refused + 1))
done
accepted=0
for args in '--meth=emseq' '--meth=em-seq' '--meth=bisulfite' '--meth=taps' '--meth=5base' '--meth=5-base' \
    '--meth=taps --meth-context=cpg' '--meth=taps --meth-scoring=genomic' \
    '--meth --meth-tags=XR,XG' '--meth --meth-tags=^XM' '--meth --meth-tags=-XM' \
    '--meth --set-as-failed=f' '--meth --set-as-failed=r'; do
    read -r -a options <<< "$args"
    for command in mem aln; do
        # /dev/null passes MEM's input check; accepted options must reach the index probe.
        rc=0
        "$bin" "$command" "${options[@]}" /nonexistent /dev/null > "$tmp/accepted.out" 2> "$tmp/accepted.err" || rc=$?
        [ "$rc" = 1 ] && [ ! -s "$tmp/accepted.out" ] || {
            printf 'FAIL: accepted-value probe exit/output\n' >&2
            exit 1
        }
        if [ "$command" = aln ]; then diagnostic='original reference index file'; else diagnostic='seed index'; fi
        grep -F "$diagnostic" "$tmp/accepted.err" > /dev/null || {
            printf 'FAIL: %s accepted values did not reach index probe: %s\n' "$command" "$args" >&2
            exit 1
        }
        grep -F 'bwa-mem3 index --meth' "$tmp/accepted.err" > /dev/null || exit 1
        accepted=$((accepted + 1))
    done
done
printf 'PASS: live shared help (%s lines), %s refusals, %s accepted value sets\n' \
    "$(wc -l < "$tmp/mem.block")" "$refused" "$accepted"
