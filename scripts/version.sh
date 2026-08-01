#!/usr/bin/env bash
# scripts/version.sh — print the canonical version string for this build.
#
# Single source of truth: version.txt at the repo root.
#
# When git is available, an informational dev suffix is appended:
#   <base>              — tarball / shallow-clone / clean-at-tag
#   <base>-dirty        — clean-at-tag with uncommitted changes
#   <base>-<sha>        — not at the tag, clean
#   <base>-<sha>-dirty  — not at the tag, with uncommitted changes
#
# "At the tag" means HEAD is exactly the commit pointed to by tag
# 'v<base>' or '<base>'. Manifest drift (HEAD at a different tag than
# version.txt indicates) is treated as "not at tag", so the resulting
# `-<sha>` suffix surfaces the drift visibly.

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"

# version.txt is the single source of truth, so a checkout that cannot supply
# it is broken rather than merely unlabelled. The old `|| echo unknown` turned
# that into a successful run printing a plausible-looking version string, which
# the Makefile then baked into the binary; fail loudly instead.
version_file="${here}/../version.txt"
if [ ! -r "$version_file" ] || [ ! -s "$version_file" ]; then
    echo "error: $version_file is missing, unreadable, or empty" >&2
    exit 1
fi
base="$(cat "$version_file")"

sha=""
dirty=""
at_target_tag=false

if command -v git > /dev/null 2>&1 \
    && git -C "${here}/.." rev-parse --git-dir > /dev/null 2>&1; then
    if sha="$(git -C "${here}/.." rev-parse --short=7 HEAD 2> /dev/null)"; then
        if ! git -C "${here}/.." diff --quiet HEAD 2> /dev/null; then
            dirty="-dirty"
        fi
        exact="$(git -C "${here}/.." describe --exact-match --tags HEAD 2> /dev/null || true)"
        case "$exact" in
            "v$base" | "$base") at_target_tag=true ;;
        esac
    else
        sha=""
    fi
fi

if [ "$at_target_tag" = "true" ] || [ -z "$sha" ]; then
    printf '%s%s\n' "$base" "$dirty"
else
    printf '%s-%s%s\n' "$base" "$sha" "$dirty"
fi
