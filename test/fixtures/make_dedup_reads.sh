#!/usr/bin/env bash
# Deterministic PE FASTQ pairs sliced from phix.fa.
#   make_dedup_reads.sh <phix.fa> <outprefix> <n_windows> <repeat>
# Emits n_windows distinct concordant pairs (97 bp, insert 300), each repeated
# <repeat> times under distinct read names. repeat=8 -> dup-rich extension jobs;
# repeat=1 -> low-dup. Byte-identical across awk implementations (no rand()).
# Generated in a single awk pass (no per-record process forks or file reopens).
#
# Each read carries a deterministic 3 bp deletion (a 100 bp window with bases
# 45-47 removed). The gap forces GAPPED extension, so the reads actually reach
# the banded Smith-Waterman kernel (and thus the dedup path under test) instead
# of resolving via the exact-match/ungapped fast path -- an exact-substring
# fixture produces zero banded-SW jobs and would make the dedup byte-identity
# test vacuous. The deletion is identical across a window's repeats, so repeats
# stay byte-identical duplicate jobs.
set -euo pipefail
FA=$1
OUT=$2
NW=$3
REP=$4
AWK="${AWK:-$(command -v mawk || command -v gawk || command -v awk)}"
# One awk pass: the main rule concatenates the FASTA sequence lines, then
# END emits every record to both mate files.
# shellcheck disable=SC2016 # $0 is awk's record, not a shell expansion
"$AWK" -v out="$OUT" -v nw="$NW" -v rep="$REP" '
    /^>/ { next }
    { seq = seq $0 }
    # Reverse complement; only ACGTacgt are mapped, anything else passes through.
    function revcomp(s,    n, j, c, r) {
        n = length(s)
        r = ""
        for (j = n; j >= 1; j--) {
            c = substr(s, j, 1)
            r = r ((c in comp) ? comp[c] : c)
        }
        return r
    }
    # Drop 3 bp (0-based chars 44..46) to force gapped extension.
    function del3(s) { return substr(s, 1, 44) substr(s, 48) }
    END {
        comp["A"] = "T"; comp["C"] = "G"; comp["G"] = "C"; comp["T"] = "A"
        comp["a"] = "t"; comp["c"] = "g"; comp["g"] = "c"; comp["t"] = "a"
        qual = sprintf("%97s", "")
        gsub(/ /, "I", qual)
        out1 = out "_1.fq"
        out2 = out "_2.fq"
        printf "" > out1
        printf "" > out2
        for (i = 0; i < nw; i++) {
            start = (i * 37) % 4900 # stride 37, stay inside 5386-300
            fwd = del3(substr(seq, start + 1, 100)) # 3 bp deletion -> banded SW
            rev = del3(revcomp(substr(seq, start + 201, 100)))
            for (r = 0; r < rep; r++) {
                printf "@dd_%05d_%02d/1\n%s\n+\n%s\n", i, r, fwd, qual > out1
                printf "@dd_%05d_%02d/2\n%s\n+\n%s\n", i, r, rev, qual > out2
            }
        }
    }
' "$FA"
