#!/usr/bin/env bash
# Deterministic PE FASTQ pairs sliced from phix.fa.
#   make_dedup_reads.sh <phix.fa> <outprefix> <n_windows> <repeat>
# Emits n_windows distinct concordant pairs (2x100, insert 300), each repeated
# <repeat> times under distinct read names. repeat=8 -> dup-rich extension jobs;
# repeat=1 -> low-dup. Byte-identical across awk implementations (no rand()).
set -euo pipefail
FA=$1; OUT=$2; NW=$3; REP=$4
SEQ=$(awk '/^>/{next}{printf "%s",$0}' "$FA")
: > "${OUT}_1.fq"; : > "${OUT}_2.fq"
for ((i=0; i<NW; i++)); do
    start=$(( (i * 37) % 4900 ))                     # stride 37, stay inside 5386-300
    fwd=${SEQ:$start:100}
    m2=${SEQ:$((start+200)):100}
    rev=$(printf '%s' "$m2" | rev | tr ACGTacgt TGCAtgca)
    for ((r=0; r<REP; r++)); do
        printf '@dd_%05d_%02d/1\n%s\n+\n%s\n' "$i" "$r" "$fwd" \
            "$(printf 'I%.0s' $(seq 100))" >> "${OUT}_1.fq"
        printf '@dd_%05d_%02d/2\n%s\n+\n%s\n' "$i" "$r" "$rev" \
            "$(printf 'I%.0s' $(seq 100))" >> "${OUT}_2.fq"
    done
done
