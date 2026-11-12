#!/usr/bin/env bash
# Generated-fixture coverage for prefix selection, minimal indexes, and safe replacement.
# Environment: BWA_MEM3 is required (candidate executable path or PATH command).
# BWA_MEM3_BASELINE optionally names a baseline executable for cross-revision parity.
# Requires python3; mawk mapping assertions and samtools BAM validation are optional.
set -Eeuo pipefail
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
trap 'fail "unexpected failure at line $LINENO"' ERR
# Normalize before entering the fixture directory, including bare PATH commands.
resolve_binary() {
    local path="$1"
    [[ $path == */* ]] || path=$(command -v "$path") || return 1
    [[ -f $path && -x $path ]] || return 1
    printf '%s/%s\n' "$(cd "$(dirname "$path")" && pwd)" "$(basename "$path")"
}
[[ -n ${BWA_MEM3:-} ]] || fail 'BWA_MEM3 must be set'
BWA_MEM3=$(resolve_binary "$BWA_MEM3") || fail 'BWA_MEM3 is not executable'
if [[ -n ${BWA_MEM3_BASELINE:-} ]]; then
    BWA_MEM3_BASELINE=$(resolve_binary "$BWA_MEM3_BASELINE") || fail 'BWA_MEM3_BASELINE is not executable'
fi
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=test/framework/mem_cli.sh
source "$HERE/../framework/mem_cli.sh"
# Reuse the record comparison without mem_cli_init's below-floor early SKIP:
# indexing checks must run even when no alignment tier is available.
WORK=$(mktemp -d) || fail 'cannot create fixture directory'
MEM_CLI_WORK=$WORK
MEM_CLI_CHECKS=0
MEM_CLI_FAILS=0
trap 'chmod -R u+w "$WORK"; rm -rf "$WORK"' EXIT
cd "$WORK"
index() { "$BWA_MEM3" index -t 1 "$@" > build.log 2>&1 || {
    cat build.log >&2
    fail "index $*"
}; }
reject() {
    if "$BWA_MEM3" index -t 1 "$@" > reject.log 2>&1; then fail "accepted $*"; fi
}
snapshot() { cksum "$1".* > "$2"; }
mkdir 'source dir' 'output dir' scratch
# A deterministic, non-repetitive sequence and proper inward-facing pairs.
python3 - << 'PY'
import random
from pathlib import Path
rng = random.Random(791)
seq = ''.join(rng.choice('ACGT') for _ in range(6000))
Path('source dir/ref.fa').write_text('>chr1\n' + seq + '\n')
rc = str.maketrans('ACGT', 'TGCA')
for chemistry in ('em', 'bs', 'taps', 'plain'):
    with open(chemistry + '.1.fq', 'w') as r1, open(chemistry + '.2.fq', 'w') as r2:
        for i, pos in enumerate(range(300, 4000, 100)):
            a, b = seq[pos:pos+100], seq[pos+200:pos+300].translate(rc)[::-1]
            if chemistry != 'plain':
                a, b = a.replace('C', 'T'), b.replace('G', 'A')
            for out, bases in ((r1, a), (r2, b)):
                out.write(f'@read{i}\n{bases}\n+\n' + 'I' * 100 + '\n')
PY
src='source dir/ref.fa'
full='output dir/full'
minimal='output dir/minimal'
reject --meth-only "$src"
grep -q -- '--meth' reject.log || fail 'missing --meth guidance'
reject --keep-meth-fasta "$src"
reject --force "$src"
# Fixed writer temporary names are unrelated inputs, never output artifacts.
for stem in temp-original temp-seed; do
    suffix=bwt.2bit.64.tmp
    [[ $stem == temp-seed ]] && suffix=meth.bwt.2bit.64.tmp
    for kind in file symlink hardlink source; do
        prefix="$stem-$kind"
        cp "$src" saved-source
        cp "$src" "$prefix.victim"
        case "$kind" in
            file) cp "$src" "$prefix.$suffix" ;;
            symlink) ln -s "$PWD/$prefix.victim" "$prefix.$suffix" ;;
            hardlink) ln "$prefix.victim" "$prefix.$suffix" ;;
            source) cp "$src" "$prefix.$suffix" ;;
        esac
        input="$src"
        [[ $kind == source ]] && input="$prefix.$suffix"
        # A secure unique-temp writer may build successfully; refusal is also safe.
        "$BWA_MEM3" index --meth -t 1 -p "$prefix" "$input" > temp.log 2>&1 || true
        "$BWA_MEM3" index --meth --force -t 1 -p "$prefix" "$input" > temp.log 2>&1 || true
        cmp saved-source "$prefix.$suffix" || fail "writer clobbered $stem $kind"
        cmp saved-source "$prefix.victim" || fail "writer clobbered victim $stem $kind"
    done
done
index --meth -p "$full" "$src"
for suffix in amb ann pac bwt.2bit.64 meth.amb meth.ann meth.pac meth.bwt.2bit.64; do
    [[ -s $full.$suffix ]] || fail "missing full $suffix"
done
[[ ! -e $full.meth.fa ]] || fail 'conversion retained by default'
# Optional cross-revision control for local validation; CI needs only the candidate.
if [[ -n ${BWA_MEM3_BASELINE:-} ]]; then
    cp "$src" baseline.fa
    "$BWA_MEM3_BASELINE" index --meth -t 1 baseline.fa > baseline.log 2>&1 || fail 'baseline index'
    # The shared FM writer must also preserve ordinary indexing bytes and umask.
    (
        umask 077
        "$BWA_MEM3" index -t 1 -p ordinary-candidate "$src"
    ) > ordinary-build.log 2>&1
    "$BWA_MEM3_BASELINE" index -t 1 -p ordinary-baseline "$src" > ordinary-baseline.log 2>&1
    for suffix in amb ann pac bwt.2bit.64; do
        cmp "ordinary-candidate.$suffix" "ordinary-baseline.$suffix" || fail "ordinary index differs: $suffix"
    done
    python3 - << 'PYPERM'
from pathlib import Path
assert Path('ordinary-candidate.bwt.2bit.64').stat().st_mode & 0o777 == 0o600
PYPERM
    for suffix in amb ann pac bwt.2bit.64 meth.amb meth.ann meth.pac meth.bwt.2bit.64; do
        cmp "$full.$suffix" "baseline.fa.$suffix" || fail "baseline index differs: $suffix"
    done
fi
index --meth --meth-only --tmp-dir scratch -p "$minimal" "$src"
grep -q 'writing fresh seed FASTA scratch/bwa-meth-' build.log || fail 'initial scratch selector ignored'
[[ $(find 'output dir' -name 'minimal.*' -type f | wc -l) -eq 7 ]] || fail 'minimal layout is not seven files'
[[ ! -e $minimal.bwt.2bit.64 && ! -e $minimal.0123 ]] || fail 'minimal contains ordinary FM-index'
[[ -z $(ls -A scratch) ]] || fail 'scratch leaked'
[[ $(find 'source dir' -type f | wc -l) -eq 1 ]] || fail 'wrote beside source'
snapshot "$full" before
reject --meth -p "$full" "$src"
snapshot "$full" after
cmp before after || fail 'rerun without force changed output'
reject --meth --meth-only --emit-unpacked-ref -p incompatible "$src"
reject --meth -p ambiguous.meth "$src"
reject --meth --force --max-memory 1K -p "$full" "$src"
grep -q 'seed index dominates' reject.log || fail 'missing seed budget diagnostic'
snapshot "$full" after
cmp before after || fail 'budget refusal changed output'
for chemistry in em bs taps; do
    args=(--meth)
    [[ $chemistry == bs ]] && args=(--meth=bisulfite)
    [[ $chemistry == taps ]] && args=(--meth=taps)
    for layout in se pe; do
        reads=("$chemistry.1.fq")
        [[ $layout == pe ]] && reads+=("$chemistry.2.fq")
        for prefix in "$full" "$minimal"; do
            "$BWA_MEM3" mem -t 1 "${args[@]}" "$prefix" "${reads[@]}" > "$prefix.sam" 2> align.log || fail "align $chemistry $layout"
            if command -v mawk > /dev/null 2>&1; then
                mawk '$1 !~ /^@/ && $3 != "*" {n++} END {exit n == 0}' "$prefix.sam" || fail 'no mappings'
            else
                echo "SKIP: mawk not on PATH; mapping assertion for $chemistry $layout $prefix not run"
            fi
        done
        mem_cli_same_file "full/minimal $chemistry $layout" "$full.sam" "$minimal.sam"
        if [[ -n ${BWA_MEM3_BASELINE:-} ]]; then
            "$BWA_MEM3_BASELINE" mem -t 1 "${args[@]}" baseline.fa "${reads[@]}" > baseline.sam 2> baseline.log || fail 'baseline align'
            mem_cli_same_file "baseline/candidate $chemistry $layout" "$full.sam" baseline.sam
        fi
    done
done
"$BWA_MEM3" mem -t 1 "$full" plain.1.fq > ordinary.sam 2> align.log || fail 'full ordinary alignment'
if command -v mawk > /dev/null 2>&1; then
    mawk '$1 !~ /^@/ && $3 != "*" {n++} END {exit n == 0}' ordinary.sam || fail 'no ordinary mappings'
else
    echo 'SKIP: mawk not on PATH; ordinary mapping assertion not run'
fi
if "$BWA_MEM3" mem -t 1 "$minimal" plain.1.fq > /dev/null 2> align.log; then fail 'ordinary mem accepted minimal'; fi
grep -q -- '--meth' align.log || fail 'minimal error lacks guidance'
printf 'stale\n' > "$full.0123"
printf 'stale\n' > "$full.meth.0123"
index --meth --force --meth-only --keep-meth-fasta -p "$full" "$src"
[[ ! -e $full.0123 && ! -e $full.meth.0123 ]] || fail 'minimal retained stale unpacked reference'
[[ ! -e $full.bwt.2bit.64 && -s $full.meth.fa ]] || fail 'full to minimal / retention'
# Stale converted FASTA must never be reused.
printf '>wrong\nAAAA\n' > "$full.meth.fa"
index --meth --force --keep-meth-fasta -p "$full" "$src"
grep -q '^>rchr1' "$full.meth.fa" || fail 'reused stale FASTA'
[[ -s $full.bwt.2bit.64 ]] || fail 'minimal to full'
printf 'stale\n' > "$full.0123"
printf 'stale\n' > "$full.meth.0123"
index --meth --force -p "$full" "$src"
[[ ! -e $full.0123 && ! -e $full.meth.0123 ]] || fail 'default retained stale unpacked reference'
[[ ! -e $full.meth.fa ]] || fail 'stale retained FASTA survived'
# Source and unrelated files must survive even with force.
printf 'unrelated\n' > unrelated
ln -s "$PWD/unrelated" collision.pac
reject --meth --force -p collision "$src"
grep -q '^unrelated$' unrelated || fail 'followed output symlink'
mkdir directory.pac
reject --meth --force -p directory "$src"
cp "$src" alias.pac
reject --meth --force -p alias alias.pac
cmp "$src" alias.pac || fail 'destroyed source alias'
cp "$src" hard-before
ln "$src" hard.pac
reject --meth --force -p hard "$src"
cmp hard-before hard.pac || fail 'destroyed hardlink'
reject --meth --tmp-dir missing -p badtmp "$src"
[[ ! -e badtmp.pac ]] || fail 'invalid scratch wrote output'
# Force a writer failure after conversion, exercising lower-level exit() cleanup.
if (
    trap '' XFSZ
    ulimit -f 32
    "$BWA_MEM3" index --meth --tmp-dir scratch -t 1 -p limited "$src"
) > limited.log 2>&1; then
    fail 'file-size limit did not fail build'
fi
grep -q 'writing fresh seed FASTA scratch/bwa-meth-' limited.log || fail 'failure never reached scratch conversion'
[[ -z $(ls -A scratch) ]] || fail 'scratch leaked on writer exit'
[[ -z $(find . -name 'limited.bwt.2bit.64.*' -o -name 'limited.meth.bwt.2bit.64.*') ]] || fail 'writer temporary leaked'
# Conversion failure occurs before replacement, unlike the later writer failure.
snapshot "$minimal" conversion-before
if (
    trap '' XFSZ
    ulimit -f 1
    "$BWA_MEM3" index --meth --force --tmp-dir scratch -t 1 -p "$minimal" "$src"
) > conversion-failure.log 2>&1; then
    fail 'conversion file limit accepted'
fi
grep -q 'failed reading or writing methylation conversion' conversion-failure.log || fail 'conversion did not fail as intended'
snapshot "$minimal" conversion-after
cmp conversion-before conversion-after || fail 'conversion failure changed outputs'
[[ -z $(ls -A scratch) ]] || fail 'conversion failure leaked scratch'
# A malformed compressed source must be rejected before replacing an index.
gzip -c "$src" > truncated.gz
python3 - << 'PYTRUNC'
from pathlib import Path
p = Path('truncated.gz')
p.write_bytes(p.read_bytes()[:-8])
PYTRUNC
snapshot "$minimal" truncated-before
reject --meth --force -p "$minimal" truncated.gz
snapshot "$minimal" truncated-after
cmp truncated-before truncated-after || fail 'truncated source changed outputs'
grep -q 'cannot read' reject.log || fail 'missing truncated-source diagnostic'
printf '>empty\n' > empty.fa
reject --meth --tmp-dir scratch -p empty empty.fa
[[ -z $(ls -A scratch) ]] || fail 'scratch leaked on failed build'
gzip -c "$src" > 'source dir/ref.fa.gz'
chmod a-w 'source dir' 'source dir/'*
index --meth --meth-only -p 'output dir/compressed' 'source dir/ref.fa.gz'
cmp "$minimal.pac" 'output dir/compressed.pac' || fail 'gzip packing differs'
chmod u+w 'source dir' 'source dir/'*
# Source .hdr wins over basename .dict, and only original contigs get metadata.
printf 'chr1\n' > "$src.alt"
printf '@SQ\tSN:chr1\tLN:6000\tAS:header-wins\n' > "$src.hdr"
printf '@SQ\tSN:chr1\tLN:6000\tAS:dictionary\n' > 'source dir/ref.dict'
index --meth --meth-only -p 'output dir/meta.fa' "$src"
cmp "$src.alt" 'output dir/meta.fa.alt' || fail 'ALT not relocated'
cmp "$src.hdr" 'output dir/meta.fa.hdr' || fail 'header precedence changed'
[[ ! -e 'output dir/meta.fa.meth.alt' && ! -e 'output dir/meta.fa.meth.hdr' ]] || fail 'seed metadata emitted'
# Compare actual relocated alignment metadata against an in-place index.
index --meth --meth-only "$src"
for prefix in "$src" 'output dir/meta.fa'; do
    "$BWA_MEM3" mem --meth -t 1 "$prefix" em.1.fq > "$prefix.sam" 2> align.log
    grep -q 'AS:header-wins' "$prefix.sam" || fail 'relocated SAM header lost metadata'
    if command -v mawk > /dev/null 2>&1; then
        mawk '$1 !~ /^@/ && $3 != "*" {n++} END {exit n == 0}' "$prefix.sam" || fail 'metadata fixture unmapped'
    else
        echo "SKIP: mawk not on PATH; metadata mapping assertion for $prefix not run"
    fi
done
mem_cli_same_file 'relocated ALT/header alignment' "$src.sam" 'output dir/meta.fa.sam'
# Without a header sidecar, ALT status is emitted from the relocated .alt.
cp "$src" alt-only.fa
printf 'chr1\n' > alt-only.fa.alt
index --meth --meth-only -p 'output dir/alt-only' alt-only.fa
"$BWA_MEM3" mem --meth -t 1 'output dir/alt-only' em.1.fq > alt-only.sam 2> align.log
grep -q 'AH:' alt-only.sam || fail 'relocated ALT annotation absent from SAM'
rm "$src.hdr"
index --meth --meth-only -p 'output dir/dict' 'source dir/ref.fa.gz'
cmp 'source dir/ref.dict' 'output dir/dict.hdr' || fail 'gzip dictionary lookup changed'
# In-place dictionary fallback must not be promoted into a stale .hdr.
cp "$src" inplace.fa
cp 'source dir/ref.dict' inplace.dict
index --meth inplace.fa
[[ ! -e inplace.fa.hdr ]] || fail 'in-place dictionary promoted to header'
printf 'A\n' >> inplace.fa
printf '@SQ\tSN:chr1\tLN:6001\tAS:updated-dictionary\n' > inplace.dict
index --meth --force inplace.fa
"$BWA_MEM3" mem --meth -t 1 inplace.fa em.1.fq > inplace.sam 2> align.log
grep -q 'AS:updated-dictionary' inplace.sam || fail 'dictionary update shadowed'
if command -v samtools > /dev/null 2>&1; then
    "$BWA_MEM3" mem --meth --bam -t 1 inplace.fa em.1.fq > inplace.bam 2> align.log
    samtools view -h inplace.bam > inplace-bam.sam
    mem_cli_same_file 'updated dictionary SAM/BAM' inplace.sam inplace-bam.sam
else
    echo 'SKIP: samtools not on PATH; dictionary-update BAM check not run'
fi
# Every conflicting metadata name is preserved and refused before output.
for suffix in dict fa.dict fa.meth.hdr fa.meth.alt; do
    printf 'unrelated metadata\n' > "conflict.$suffix"
    cp "conflict.$suffix" saved-metadata
    for forced in '' --force; do
        reject --meth ${forced:+"$forced"} -p conflict.fa "$src"
        grep -q 'conflicting metadata' reject.log || fail 'missing metadata diagnostic'
        cmp saved-metadata "conflict.$suffix" || fail 'metadata clobbered'
        [[ ! -e conflict.fa.pac ]] || fail 'metadata refusal wrote output'
    done
    rm "conflict.$suffix"
done
# Existing unreadable metadata must never be mistaken for absence.
if [[ $(id -u) != 0 ]]; then
    cp "$src" unreadable.fa
    for suffix in alt hdr dict; do
        path="unreadable.fa.$suffix"
        [[ $suffix == dict ]] && path=unreadable.dict
        printf 'metadata\n' > "$path"
        chmod 000 "$path"
        snapshot "$minimal" unreadable-before
        reject --meth --force -p "$minimal" unreadable.fa
        snapshot "$minimal" unreadable-after
        cmp unreadable-before unreadable-after || fail 'unreadable sidecar changed outputs'
        grep -q 'metadata' reject.log || fail 'missing unreadable metadata diagnostic'
        chmod u+rw "$path"
        rm "$path"
    done
else
    echo 'SKIP: unreadable sidecar permission test requires non-root user'
fi
# Explicit tmp-dir wins over both environment variables; BWA_INDEX_TMPDIR wins over TMPDIR.
BWA_INDEX_TMPDIR=missing TMPDIR=missing index --meth --meth-only --tmp-dir scratch -p tmpcli "$src"
grep -q 'writing fresh seed FASTA scratch/bwa-meth-' build.log || fail 'CLI scratch precedence ignored'
BWA_INDEX_TMPDIR=scratch TMPDIR=missing index --meth --meth-only -p tmpenv "$src"
grep -q 'writing fresh seed FASTA scratch/bwa-meth-' build.log || fail 'BWA_INDEX_TMPDIR precedence ignored'
(
    unset BWA_INDEX_TMPDIR
    TMPDIR=scratch index --meth --meth-only -p tmpsystem "$src"
)
grep -q 'writing fresh seed FASTA scratch/bwa-meth-' build.log || fail 'TMPDIR ignored'
[[ -z $(ls -A scratch) ]] || fail 'scratch leaked after success'
[[ $MEM_CLI_CHECKS -gt 0 && $MEM_CLI_FAILS -eq 0 ]] || fail 'SAM comparison failed'
echo 'PASS: meth_index_output (layouts, alignment parity, safety, metadata, scratch, gzip)'
