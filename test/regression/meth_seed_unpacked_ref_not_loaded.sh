#!/usr/bin/env bash
# test/regression/meth_seed_unpacked_ref_not_loaded.sh
#
# Regression (D3 memory): the SEED index's unpacked reference `<ref>.meth.0123`
# is neither BUILT nor LOADED. D3 seeds in the doubled `.meth` FM-index but
# extends/scores against the ORIGINAL reference (meth_orig_*), so the seed
# `.0123` is dead weight (~13 GB on hg38, never read). `index --meth` therefore
# does not emit it, and `mem --meth` does not load it.
#
# Asserts:
#   1. `index --meth` builds the seed FM-index/bns/pac but NOT `<ref>.meth.0123`.
#   2. The ORIGINAL `<ref>.0123` (the real extension target) IS built.
#   3. `mem --meth` runs to a valid, non-empty BAM without the seed `.0123`, and
#      the output is deterministic.
#
# RED on a binary that builds + loads the seed `.0123`: assertion 1 fails (the
# file is present). GREEN once the seed `.0123` is dropped from build + load.
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
set -euo pipefail
: "${BWA_MEM3:?BWA_MEM3 must be set}"
command -v samtools >/dev/null 2>&1 || { echo "SKIP: samtools not on PATH (--meth emits BAM)"; exit 0; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT; cd "$WORK"
fail() { echo "FAIL: $*" >&2; exit 1; }
rc() { printf '%s' "$1" | rev | tr ACGTacgt TGCAtgca; }

# Deterministic 1500 bp reference (same generator/seed as meth_output_integrity).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG

printf '>chrA\n%s\n' "$REF" > ref.fa
Q=$(printf 'I%.0s' $(seq 1 60))
"$BWA_MEM3" index --meth ref.fa >/dev/null 2>&1 || fail "index --meth nonzero exit"

# 1. The seed unpacked reference must NOT be built; the seed FM-index/bns/pac must.
[ ! -e ref.fa.meth.0123 ] || fail "seed ref.fa.meth.0123 was built; it must not be (never read in --meth)"
for f in ref.fa.meth.bwt.2bit.64 ref.fa.meth.ann ref.fa.meth.amb ref.fa.meth.pac; do
  [ -s "$f" ] || fail "expected seed index file $f to be built"
done
# 2. The ORIGINAL unpacked reference (the real extension target) must be built.
[ -s ref.fa.0123 ] || fail "expected original ref.fa.0123 (the --meth extension target) to be built"

# Proper FR pairs derived from the reference (exact windows so they place cleanly).
: > r1.fq; : > r2.fq
mk_pair() { local r1 r2; r1=${REF:$2:60}; r2=$(rc "${REF:$3:60}")
  printf '@%s\n%s\n+\n%s\n' "$1" "$r1" "$Q" >> r1.fq
  printf '@%s\n%s\n+\n%s\n' "$1" "$r2" "$Q" >> r2.fq; }
mk_pair p1 100 300; mk_pair p2 500 720; mk_pair p3 900 1140

run_meth() { "$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa r1.fq r2.fq > "$1" 2>/dev/null; }

# 3. mem --meth runs without the seed `.0123`, emits a valid non-empty BAM,
#    and is deterministic.
run_meth a.bam || fail "mem --meth nonzero exit (must not need the seed .0123)"
samtools quickcheck a.bam || fail "a.bam invalid"
samtools view a.bam > a.records || fail "samtools view a.bam failed"
[ -s a.records ] || fail "mem --meth produced no alignment records"
run_meth b.bam || fail "mem --meth nonzero exit on second run"
samtools quickcheck b.bam || fail "b.bam invalid"
samtools view b.bam > b.records || fail "samtools view b.bam failed"
diff -q a.records b.records >/dev/null 2>&1 || fail "mem --meth output is non-deterministic"

echo "PASS: meth_seed_unpacked_ref_not_loaded (seed .meth.0123 is neither built nor loaded; mem --meth works without it)"
