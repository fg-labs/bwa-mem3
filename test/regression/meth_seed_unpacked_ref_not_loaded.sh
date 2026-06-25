#!/usr/bin/env bash
# test/regression/meth_seed_unpacked_ref_not_loaded.sh
#
# Regression (D3 memory): the SEED index's unpacked reference (`<ref>.meth.0123`)
# is NOT consumed by `mem --meth`. D3 seeds in the doubled `.meth` FM-index but
# extends/scores against the ORIGINAL reference (meth_orig_*), so the seed
# `.0123` is dead weight (~13 GB on hg38). This test pins the contract that it is
# never loaded, by deleting it after indexing and asserting `mem --meth` still
# runs and produces a BYTE-IDENTICAL BAM.
#
# RED on a binary that still loads the seed `.0123`: the second `mem --meth`
# exits non-zero (load_ref_string fails on the missing file).
# GREEN once the seed `.0123` load is gated off in --meth: identical output.
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
[ -f ref.fa.meth.0123 ] || fail "expected seed index ref.fa.meth.0123 to exist after indexing"

# Proper FR pairs derived from the reference (exact windows so they place
# cleanly and deterministically; exercises seeding + original-ref extension).
mk_pair() {  # $1=name $2=fwd-start(0-based) $3=rev-start
  local r1 r2
  r1=${REF:$2:60}
  r2=$(rc "${REF:$3:60}")
  printf '@%s\n%s\n+\n%s\n' "$1" "$r1" "$Q" >> r1.fq
  printf '@%s\n%s\n+\n%s\n' "$1" "$r2" "$Q" >> r2.fq
}
: > r1.fq; : > r2.fq
mk_pair p1 100 300
mk_pair p2 500 720
mk_pair p3 900 1140

run_meth() { # $1=out.bam
  "$BWA_MEM3" mem --meth --meth-scoring genomic -t 1 ref.fa r1.fq r2.fq > "$1" 2>/dev/null
}

# 1. Baseline run with the full index present (the golden output).
run_meth with.bam || fail "mem --meth nonzero exit with full index"
samtools quickcheck with.bam || fail "with.bam invalid"
samtools view with.bam > with.records || fail "samtools view with.bam failed"
[ -s with.records ] || fail "with.bam produced no alignment records"

# 2. Delete the SEED unpacked reference and re-run. A binary that does not
#    consume it must succeed; one that loads it will exit non-zero here.
rm -f ref.fa.meth.0123
run_meth without.bam || fail "mem --meth exited non-zero after removing ref.fa.meth.0123 (the seed .0123 must not be loaded in --meth)"
samtools quickcheck without.bam || fail "without.bam invalid"
samtools view without.bam > without.records || fail "samtools view without.bam failed"

# 3. Output must be byte-identical: removing a never-read buffer cannot change a thing.
diff -q with.records without.records >/dev/null 2>&1 \
  || fail "BAM records changed after removing the seed .0123 (output must be identical)"

echo "PASS: meth_seed_unpacked_ref_not_loaded (seed .meth.0123 is not consumed by mem --meth; output byte-identical without it)"
