# Unit-test fixtures

Committed inputs for the five C++ unit tests under `test/`:

| Binary | Input file | Format |
|---|---|---|
| `fmi_test` | `reads.fa` | ACGT FASTA, 10 × 50 bp |
| `smem2_test` | `smem2_input.txt` | Alternating lines: 50-char `0/1/2/3` sequence, then `<query_pos> <min_intv>` |
| `bwt_seed_strategy_test` | `bwt_seed_input.fa` | FASTA with 50-char `0/1/2/3` sequences |
| `sa2ref_test` | `sa2ref_input.txt` | Lines of `"<k>, <s>"` — suffix-array offset + count |
| `xeonbsw` | `pairs.txt` | Triplets: `<h0>\n<ref 0/1/2/3>\n<query 0/1/2/3>\n` |

All four SMEM/seed tests share `phix.fa` (phiX174, NC_001422.1, 5386 bp) as the reference. The FMI index is built on demand by the harness — not committed.

`synthetic_1mb.fa` is a 1 Mbp ACGT-only FASTA used by the libsais byte-diff, determinism, and memory-budget tests; the matching `.0123` and `.bwt.2bit.64` baselines under `baselines/` are derived from it and only stay reproducible when both source and baselines are committed together.

## Regenerating the fixtures

Fixtures are derived from `phix.fa` (committed). To regenerate after a format change:

```sh
# Rebuild phix.fa from NCBI (only if you think it has drifted — it should not)
curl -sSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta" | \
  mawk 'BEGIN{RS=">"; ORS=""} NR>1 { h=$1; sub(/[^\n]+/, "", $0); gsub(/\n/, "", $0); printf ">%s\n", h; for(i=1;i<=length($0);i+=70) printf "%s\n", substr($0, i, 70) }' > phix.fa

# reads.fa — 10 × 50 bp ACGT, stride 500
mawk 'BEGIN{RS=">"; ORS=""} NR==2 { sub(/[^\n]+/, "", $0); gsub(/\n/, "", $0); for(i=0;i<10;i++){ off = i * 500; if (off+50 > length($0)) break; printf ">r%d\n%s\n", i+1, substr($0, off+1, 50) } }' phix.fa > reads.fa

# smem2_input.txt — 5 × 50 bp 0/1/2/3-encoded, query_pos=25, min_intv=1
mawk 'BEGIN{RS=">"; ORS=""} NR==2 { sub(/[^\n]+/, "", $0); gsub(/\n/, "", $0); for(i=0;i<5;i++){ off = i * 700; if (off+50 > length($0)) break; s = substr($0, off+1, 50); gsub(/A/, "0", s); gsub(/C/, "1", s); gsub(/G/, "2", s); gsub(/T/, "3", s); printf "%s\n25 1\n", s } }' phix.fa > smem2_input.txt

# bwt_seed_input.fa — 5 × 50 bp 0/1/2/3-encoded FASTA
mawk 'BEGIN{RS=">"; ORS=""} NR==2 { sub(/[^\n]+/, "", $0); gsub(/\n/, "", $0); for(i=0;i<5;i++){ off = i * 900; if (off+50 > length($0)) break; s = substr($0, off+1, 50); gsub(/A/, "0", s); gsub(/C/, "1", s); gsub(/G/, "2", s); gsub(/T/, "3", s); printf ">r%d\n%s\n", i+1, s } }' phix.fa > bwt_seed_input.fa

# sa2ref_input.txt — fixed small k values
cat > sa2ref_input.txt <<EOF
0, 1
100, 1
500, 1
1000, 1
2000, 1
EOF

# pairs.txt — 3 hand-crafted triplets (do not regenerate unless you update test/run_unit_tests.sh)
```

## Running locally

```sh
cd ..                      # back to repo root
make                       # or `make arm64` on Apple Silicon
bash test/run_unit_tests.sh
```

The harness exits 0 when every binary runs to completion and emits non-empty output.
