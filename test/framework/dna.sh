#!/usr/bin/env bash
# test/framework/dna.sh
#
# Source-only DNA string helpers for regression scripts that build reads from a
# reference, so that every such script reverse-complements the same way. Keeps
# to bash 3.2 (macOS /bin/bash).
#
#     source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/dna.sh"
#     r2="$(dna_revcomp "$fwd")"

# dna_revcomp <seq>
#   The reverse complement of <seq>, case preserved: A/C/G/T and a/c/g/t swap
#   with their complements, and every other character (N, '-') is kept.
dna_revcomp() {
    printf '%s\n' "$1" | rev | tr 'ACGTacgt' 'TGCAtgca'
}
