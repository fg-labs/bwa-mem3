#!/usr/bin/env bash
# Ordinary indexing refuses existing output artifacts unless --force is given.
# Forced writes are not transactional; no test here claims concurrent-reader safety.
# Env: BWA_MEM3 (required executable). Fixtures are generated in private scratch.
set -euo pipefail
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
[[ -n ${BWA_MEM3:-} ]] || fail "BWA_MEM3 must be set"
case "$BWA_MEM3" in
    /*) BIN=$BWA_MEM3 ;;
    *) BIN="$PWD/$BWA_MEM3" ;;
esac
[[ -x $BIN ]] || fail "BWA_MEM3 is not executable"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
printf '>chrA\n' > a.fa
printf '>chrB\n' > b.fa
for _ in {1..100}; do
    printf 'ACGTTGCAAGCTTAGCGATCCGATGCTAGCAT\n' >> a.fa
    printf 'TGCAACGTCGATTACGCATGCCTAGGATCGTAACG\n' >> b.fa
done
index() { "$BIN" index -t 1 "$@" > build.out 2> build.err || {
    cat build.err >&2
    fail "index $*"
}; }
reject() {
    if "$BIN" index -t 1 "$@" > reject.out 2> reject.err; then fail "accepted index $*"; fi
    grep -qF 'use --force' reject.err || {
        cat reject.err >&2
        fail "missing --force guidance"
    }
}
index --emit-unpacked-ref a.fa
mkdir saved
for suffix in amb ann pac bwt.2bit.64 0123; do cp "a.fa.$suffix" "saved/$suffix"; done
cp a.fa saved/input
reject a.fa
reject -p a.fa b.fa
for suffix in amb ann pac bwt.2bit.64 0123; do cmp "a.fa.$suffix" "saved/$suffix" || fail "refusal changed $suffix"; done
cmp a.fa saved/input || fail "refusal changed the input FASTA"
echo 'PASS: same-reference and changed-reference overwrites refused without mutation'

index --force --emit-unpacked-ref a.fa
for suffix in amb ann pac bwt.2bit.64 0123; do cmp "a.fa.$suffix" "saved/$suffix" || fail "forced same-reference build changed $suffix"; done
index --force --emit-unpacked-ref -p a.fa b.fa
index --emit-unpacked-ref -p fresh b.fa
for suffix in amb ann pac bwt.2bit.64 0123; do cmp "a.fa.$suffix" "fresh.$suffix" || fail "forced replacement differs from fresh $suffix"; done
cmp a.fa saved/input || fail "forced rebuild changed input FASTA"
echo 'PASS: forced same-reference and explicit-prefix changed-reference builds match fresh output'

# Any one artifact is enough, including a partial index and an optional .0123.
for suffix in amb ann pac bwt.2bit.64 0123; do
    mkdir "partial-$suffix"
    printf 'existing artifact\n' > "partial-$suffix/ref.$suffix"
    cp "partial-$suffix/ref.$suffix" expected
    reject -p "partial-$suffix/ref" b.fa
    cmp "partial-$suffix/ref.$suffix" expected || fail "partial $suffix changed"
    [[ $(find "partial-$suffix" -type f | wc -l) -eq 1 ]] || fail "partial refusal wrote files"
    index --force -p "partial-$suffix/ref" b.fa
    for core in amb ann pac bwt.2bit.64; do cmp "partial-$suffix/ref.$core" "fresh.$core" || fail "forced partial $suffix differs"; done
done
echo 'PASS: partial indexes require --force and can then be rebuilt'

# Validate every destination before modifying any earlier artifact. Force must
# not follow symlinks, overwrite a source alias, or replace unrelated sidecars.
for kind in symlink dangling hardlink directory input; do
    mkdir "$kind"
    printf 'keep this artifact\n' > "$kind/ref.ann"
    cp "$kind/ref.ann" keep-ann
    case "$kind" in
        symlink) ln -s "$WORK/b.fa" "$kind/ref.pac" ;;
        dangling) ln -s "$WORK/missing" "$kind/ref.pac" ;;
        hardlink) ln b.fa "$kind/ref.pac" ;;
        directory) mkdir "$kind/ref.pac" ;;
        input) cp b.fa "$kind/ref.pac" ;;
    esac
    input=b.fa
    [[ $kind != input ]] || input="$kind/ref.pac"
    cp "$input" keep-input
    if "$BIN" index --force -t 1 -p "$kind/ref" "$input" > unsafe.out 2> unsafe.err; then fail "accepted unsafe $kind output"; fi
    grep -qF 'unsafe output collision' unsafe.err || {
        cat unsafe.err >&2
        fail "wrong unsafe-output diagnostic"
    }
    cmp "$kind/ref.ann" keep-ann || fail "unsafe preflight changed an earlier output"
    cmp "$input" keep-input || fail "unsafe preflight changed the input"
done
# The unsafe destination comes after safe, replaceable .pac/.ann/.amb files: a
# writer that validated each artifact only as it reached it would replace them first.
mkdir late-unsafe
for suffix in pac ann amb; do printf 'earlier %s artifact\n' "$suffix" > "late-unsafe/ref.$suffix"; done
cp -R late-unsafe keep-late
ln -s "$WORK/a.fa" late-unsafe/ref.bwt.2bit.64
cp a.fa keep-input
if "$BIN" index --force -t 1 -p late-unsafe/ref b.fa > unsafe.out 2> unsafe.err; then fail "accepted unsafe .bwt.2bit.64 output"; fi
grep -qF 'unsafe output collision' unsafe.err || {
    cat unsafe.err >&2
    fail "wrong unsafe-output diagnostic for .bwt.2bit.64"
}
for suffix in pac ann amb; do cmp "late-unsafe/ref.$suffix" "keep-late/ref.$suffix" || fail "unsafe .bwt.2bit.64 preflight changed earlier $suffix"; done
[[ $(find late-unsafe -type f | wc -l) -eq 3 ]] || fail "unsafe .bwt.2bit.64 preflight wrote files"
cmp a.fa keep-input || fail "unsafe .bwt.2bit.64 preflight changed the symlink target"
mkdir sidecars
printf 'unrelated ALT metadata\n' > sidecars/ref.alt
printf 'unrelated header metadata\n' > sidecars/ref.hdr
cp sidecars/ref.alt keep-alt
cp sidecars/ref.hdr keep-hdr
index -p sidecars/ref b.fa
index --force -p sidecars/ref b.fa
cmp sidecars/ref.alt keep-alt || fail "ordinary indexing replaced .alt"
cmp sidecars/ref.hdr keep-hdr || fail "ordinary indexing replaced .hdr"
echo 'PASS: forced indexing rejects aliases/symlinks/directories and preserves sidecars'
echo 'PASS: index_overwrite'
