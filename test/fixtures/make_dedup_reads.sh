#!/usr/bin/env bash
# Deterministic PE FASTQ pairs sliced from phix.fa.
#   make_dedup_reads.sh <phix.fa> <outprefix> <n_windows> <repeat>
# Emits n_windows distinct concordant pairs (97 bp, insert 300), each repeated
# <repeat> times under distinct read names. repeat=8 -> dup-rich extension jobs;
# repeat=1 -> low-dup. Byte-identical across awk implementations (no rand()).
#
# Each read carries a deterministic 3 bp deletion (a 100 bp window with bases
# 45-47 removed). The gap forces GAPPED extension, so the reads actually reach
# the banded Smith-Waterman kernel (and thus the dedup path under test) instead
# of resolving via the exact-match/ungapped fast path -- an exact-substring
# fixture produces zero banded-SW jobs and would make the dedup byte-identity
# test vacuous. The deletion is identical across a window's repeats, so repeats
# stay byte-identical duplicate jobs.
set -euo pipefail
FA=$1; OUT=$2; NW=$3; REP=$4
SEQ=$(awk '/^>/{next}{printf "%s",$0}' "$FA")
QUAL=$(printf 'I%.0s' $(seq 97))
: > "${OUT}_1.fq"; : > "${OUT}_2.fq"
for ((i=0; i<NW; i++)); do
    start=$(( (i * 37) % 4900 ))                     # stride 37, stay inside 5386-300
    fwd=${SEQ:$start:100};              fwd="${fwd:0:44}${fwd:47}"      # 3 bp deletion -> banded SW
    m2=${SEQ:$((start+200)):100}
    rev=$(printf '%s' "$m2" | rev | tr ACGTacgt TGCAtgca); rev="${rev:0:44}${rev:47}"
    for ((r=0; r<REP; r++)); do
        printf '@dd_%05d_%02d/1\n%s\n+\n%s\n' "$i" "$r" "$fwd" "$QUAL" >> "${OUT}_1.fq"
        printf '@dd_%05d_%02d/2\n%s\n+\n%s\n' "$i" "$r" "$rev" "$QUAL" >> "${OUT}_2.fq"
    done
done
