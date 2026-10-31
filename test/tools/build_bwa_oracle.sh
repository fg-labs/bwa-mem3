#!/bin/bash
# Build the stock bwa 0.7.19 reference ("oracle") binary used to check the
# bwa-mem3 aln port against upstream.
#
# usage: build_bwa_oracle.sh <out-dir>
#
# Prints the absolute path of the built `bwa` on stdout (all other output goes
# to stderr). Idempotent: an existing <out-dir>/bwa-0.7.19/bwa that reports
# Version 0.7.19 is reused. The tarball is verified against a pinned sha256
# before it is unpacked.

set -euo pipefail

VERSION="0.7.19"
URL="https://github.com/lh3/bwa/archive/refs/tags/v${VERSION}.tar.gz"
# sha256 of the v0.7.19 tag tarball, computed once on 2026-10-30.
SHA256="cdff5db67652c5b805a3df08c4e813a822c65791913eccfb3cf7d528588f37bc"

if (($# != 1)); then
    echo "usage: ${BASH_SOURCE[0]##*/} <out-dir>" >&2
    exit 2
fi

mkdir -p "$1"
out="$(cd "$1" && pwd)"
src="$out/bwa-$VERSION"
bwa="$src/bwa"

if [[ -x "$bwa" ]] && ("$bwa" 2>&1 || true) | grep -q "Version: $VERSION"; then
    echo "$bwa"
    exit 0
fi

sha256_of() {
    if command -v sha256sum > /dev/null 2>&1; then
        sha256sum "$1" | cut -d' ' -f1
    else
        shasum -a 256 "$1" | cut -d' ' -f1
    fi
}

tarball="$out/v$VERSION.tar.gz"
curl -sSfL -o "$tarball" "$URL" >&2
if [[ "$(sha256_of "$tarball")" != "$SHA256" ]]; then
    echo "FAIL: bwa-$VERSION tarball sha256 mismatch" >&2
    exit 1
fi

rm -rf "$src"
tar -xzf "$tarball" -C "$out"
(cd "$src" && make CC="${CC:-cc}" CFLAGS="-g -Wall -Wno-unused-function -O2 -ffp-contract=off" >&2)

if ! ("$bwa" 2>&1 || true) | grep -q "Version: $VERSION"; then
    echo "FAIL: built bwa does not report Version: $VERSION" >&2
    exit 1
fi
echo "$bwa"
