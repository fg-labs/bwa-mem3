#!/usr/bin/env bash
# test/regression/meth_seed_index.sh
#
# Regression: D3 `index --meth` builds a dual index — the ORIGINAL-alphabet
# index at the bare prefix plus a converted SEED FM-index at `<ref>.meth.*`
# with f/r-doubled contigs (rchrX = G->A, fchrX = C->T) — and a strand-projected
# read seeds to the correct f/r seed contig at the correct position. That locks
# in the seed->original remap arithmetic used by `mem --meth`:
#
#   orig_tid   = seed_rid / 2     (contig order r0,f0,r1,f1,...)
#   hypothesis = seed_rid & 1     (odd = f = OT [C->T];  even = r = OB [G->A])
#   orig_pos   = seed local pos   (f/r contigs are forward projections)
#
# The seed index is itself a normal bwa index, so we validate by aligning a
# projected read against it with plain `mem` and asserting the f/r contig + pos.
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# Two contigs with mixed C/T/G content so projection is non-trivial.
CHRA=ACGTACGTCCGGTTAACCGGATCGATCGCATGCATGCCGTACGTACGTACGTAACCGGTT
CHRB=TTGGCCAATTGGCCAAGCGCGCGATATATCGCGCGCGTTAACCGGATCGATGCATGCATG
printf '>chrA\n%s\n>chrB\n%s\n' "$CHRA" "$CHRB" > ref.fa

fail() { echo "FAIL: $*" >&2; exit 1; }

# --- A2: dual index emitted ------------------------------------------------
"$BWA_MEM3" index --meth ref.fa >/dev/null 2>&1 || fail "index --meth nonzero exit"
for f in ref.fa.pac ref.fa.bwt.2bit.64 ref.fa.meth.pac ref.fa.meth.bwt.2bit.64 ref.fa.meth.ann; do
    [ -s "$f" ] || fail "missing/empty index file: $f"
done

# Seed index has the f/r-doubled contigs in r0,f0,r1,f1 order.
# The bwa .ann format is two lines per contig after a one-line header:
#   "<gi> <name> <anno>"  then  "<offset> <len> <n_ambs>"
# so contig NAMES land on even line numbers (NR 2,4,6,...) in field $2 — hence
# the (NR%2)==0 filter (NOT one record per line).
ann_names=$(mawk 'NR>1 && (NR%2)==0 {print $2}' ref.fa.meth.ann | tr '\n' ' ')
[ "$ann_names" = "rchrA fchrA rchrB fchrB " ] \
    || fail "seed .meth.ann contig order: got '$ann_names' (want 'rchrA fchrA rchrB fchrB ')"

# --- A2: original index is byte-identical to a plain index ------------------
cp ref.fa refp.fa
"$BWA_MEM3" index refp.fa >/dev/null 2>&1 || fail "plain index nonzero exit"
for ext in pac ann amb bwt.2bit.64; do
    cmp -s "ref.fa.$ext" "refp.fa.$ext" || fail "original --meth index differs from plain index (.$ext)"
done

# --- remap validation: projected read -> correct f/r contig + pos ----------
# The `.meth` seed index has no `.0123` (never read by `mem --meth`), so plain
# `mem` cannot be pointed straight at the seed prefix. Probe the seed FM content
# via a plain index of the converted seed FASTA `ref.fa.meth.fa` (byte-identical
# doubled f/r text), which `index --meth` left on disk. Plain `mem` pac-fetches
# its bases from `.pac`, so no `.0123` is needed. Same contigs/positions as the
# seed `ref.fa.meth.*`.
"$BWA_MEM3" index ref.fa.meth.fa >/dev/null 2>&1 || fail "probe index of ref.fa.meth.fa nonzero exit"
# helper: align one read (seq) against the seed FM content, assert RNAME/POS.
check() {  # $1=label  $2=seq  $3=want_rname  $4=want_pos
    printf '@r\n%s\n+\n%s\n' "$2" "$(printf 'I%.0s' $(seq 1 ${#2}))" > q.fq
    line=$("$BWA_MEM3" mem ref.fa.meth.fa q.fq 2>/dev/null | mawk '$1!~/^@/{print; exit}')
    rn=$(echo "$line" | cut -f3); pos=$(echo "$line" | cut -f4)
    [ "$rn" = "$3" ] && [ "$pos" = "$4" ] \
        || fail "$1: mapped $rn:$pos, want $3:$4"
}

ct() { echo "$1" | tr C T; }   # OT projection (R1)
ga() { echo "$1" | tr G A; }   # OB projection (R2)

A40=${CHRA:0:40}; B40=${CHRB:0:40}
check "chrA OT (C->T)->fchrA" "$(ct "$A40")" fchrA 1   # rid 1 (odd) -> chrA OT
check "chrA OB (G->A)->rchrA" "$(ga "$A40")" rchrA 1   # rid 0 (even) -> chrA OB
check "chrB OT (C->T)->fchrB" "$(ct "$B40")" fchrB 1   # rid 3 (odd) -> chrB OT
check "chrB OB (G->A)->rchrB" "$(ga "$B40")" rchrB 1   # rid 2 (even) -> chrB OB

# --- A1: actionable error when the .meth seed index is absent ---------------
# refp.fa has only a plain index (no .meth.*); mem --meth must fail with guidance.
printf '@r\nACGT\n+\nIIII\n' > rr.fq
if "$BWA_MEM3" mem --meth refp.fa rr.fq >/dev/null 2>err.txt; then
    fail "mem --meth on a non-meth index should exit nonzero"
fi
grep -q 'index --meth' err.txt || fail "missing-.meth error not actionable: $(cat err.txt)"

echo "PASS: meth_seed_index (dual index + f/r layout + seed->original remap + missing-index error)"
