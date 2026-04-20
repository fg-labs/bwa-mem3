#!/usr/bin/env bash
# test/accel_parity.sh
# R1 parity check: verify that `bwa-mem2 mem --accel-cache ...` produces
# byte-identical SAM (modulo thread ordering) to stock `bwa-mem2 mem`.
#
# Expects to be run from the repository root (or test/ dir).
# Build artifacts used:
#   ../bwa-mem2 (or ./bwa-mem2 if run from repo root)
#
# Supply a reference + FASTQ or let the script synthesize one on phiX.

set -euo pipefail

# Resolve paths.
if [ -x "./bwa-mem2" ]; then
    BM2="./bwa-mem2"
elif [ -x "../bwa-mem2" ]; then
    BM2="../bwa-mem2"
else
    echo "ERROR: cannot find bwa-mem2 binary (expected ./bwa-mem2 or ../bwa-mem2)" >&2
    exit 2
fi

REF="${REF:-}"
FQ="${FQ:-}"
WORK="${WORK:-/tmp/accel_parity}"

mkdir -p "$WORK"

# If no REF given, synthesize a phiX setup under $WORK.
if [ -z "$REF" ]; then
    # Find phiX FASTA from common locations, or download from NCBI.
    PHIX=""
    for cand in \
        "$HOME/work/references/phix/phix.fasta" \
        "/tmp/phix_bm2/phix.fasta" \
        "/usr/share/phix/phix.fasta"; do
        if [ -f "$cand" ]; then PHIX="$cand"; break; fi
    done
    if [ -z "$PHIX" ]; then
        echo "ERROR: no phiX FASTA available and REF not set" >&2
        exit 2
    fi
    cp "$PHIX" "$WORK/phix.fasta"
    REF="$WORK/phix.fasta"
    # Ensure bwa-mem2 index exists.
    if [ ! -f "$REF.bwt.2bit.64" ]; then
        "$BM2" index "$REF" >/dev/null 2>&1
    fi
    # Ensure .fai exists (build-accel's BED loader needs it; samtools or
    # bwa-mem2's own indexing leaves one behind on most workflows).
    if [ ! -f "$REF.fai" ]; then
        if command -v samtools >/dev/null 2>&1; then
            samtools faidx "$REF"
        else
            # Minimal fai: single-contig.
            python3 - "$REF" <<'PY'
import sys
ref = sys.argv[1]
name = None
seq_len = 0
line_blen = None
line_len = None
offset = None
with open(ref, 'rb') as f:
    while True:
        pos = f.tell()
        line = f.readline()
        if not line: break
        if line.startswith(b'>'):
            name = line[1:].split()[0].decode()
            offset = f.tell()
            continue
        clean = line.rstrip(b'\r\n')
        if line_blen is None:
            line_blen = len(clean)
            line_len = len(line)
        seq_len += len(clean)
with open(ref + '.fai', 'w') as fo:
    fo.write(f"{name}\t{seq_len}\t{offset}\t{line_blen}\t{line_len}\n")
PY
        fi
    fi
fi

# Synthesize a small read set if no FQ supplied: 500 random 150bp reads
# sampled from the reference.
if [ -z "$FQ" ]; then
    FQ="$WORK/reads.fq"
    if [ ! -f "$FQ" ]; then
        python3 - "$REF" "$FQ" <<'PY'
import random, sys
ref_path, fq_path = sys.argv[1], sys.argv[2]
random.seed(7)
seq = []
with open(ref_path) as f:
    for line in f:
        if line.startswith('>'): continue
        seq.append(line.strip())
seq = ''.join(seq)
with open(fq_path, 'w') as out:
    for i in range(500):
        s = random.randint(0, len(seq) - 150)
        out.write(f"@r{i}\n{seq[s:s+150]}\n+\n{'I'*150}\n")
PY
    fi
fi

# BED covering the whole reference (so the accelerator cache is fully
# applicable).
if [ ! -f "$WORK/panel.bed" ]; then
    awk 'BEGIN{OFS="\t"} {print $1,0,$2}' "$REF.fai" > "$WORK/panel.bed"
fi

# Build the cache.
if [ ! -f "$WORK/panel.cache" ]; then
    "$BM2" build-accel --ref "$REF" --bed "$WORK/panel.bed" \
        --out "$WORK/panel.cache" 2>/dev/null
fi

echo "[parity] stock run"
"$BM2" mem -t 4 "$REF" "$FQ" 2>"$WORK/stock.err" > "$WORK/stock.sam"
echo "[parity] accel run"
"$BM2" mem -t 4 --accel-cache "$WORK/panel.cache" "$REF" "$FQ" \
    2>"$WORK/accel.err" > "$WORK/accel.sam"

# Confirm the accelerator actually loaded.
if ! grep -q '\[accel\] enabled' "$WORK/accel.err"; then
    echo "ERROR: accelerator did not load (see $WORK/accel.err)" >&2
    exit 3
fi

# Diff headers (strip the @PG line which embeds argv).
grep -E '^@' "$WORK/stock.sam" | grep -v '^@PG' > "$WORK/stock.hdr"
grep -E '^@' "$WORK/accel.sam" | grep -v '^@PG' > "$WORK/accel.hdr"
if ! diff -q "$WORK/stock.hdr" "$WORK/accel.hdr" >/dev/null; then
    echo "ERROR: SAM headers differ" >&2
    diff "$WORK/stock.hdr" "$WORK/accel.hdr" | head -20
    exit 4
fi

# Diff alignment records after stable sort (avoids thread-order variance).
grep -v '^@' "$WORK/stock.sam" | sort > "$WORK/stock.body"
grep -v '^@' "$WORK/accel.sam" | sort > "$WORK/accel.body"
if ! diff -q "$WORK/stock.body" "$WORK/accel.body" >/dev/null; then
    echo "ERROR: alignment records differ" >&2
    diff "$WORK/stock.body" "$WORK/accel.body" | head -40
    exit 5
fi

echo "[parity] PARITY OK — $(wc -l < "$WORK/accel.body") records match"
