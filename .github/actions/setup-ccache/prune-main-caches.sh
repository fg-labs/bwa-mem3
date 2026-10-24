#!/usr/bin/env bash
# Delete main's superseded ccache entries for one configuration.
#
# Usage: prune-main-caches.sh <key-prefix> [<saved-key>]
#
# Every push to main saves a fresh ccache entry per configuration under
# <key-prefix><sha> (see action.yml). PRs restore only the newest one, so the
# older entries are dead weight against the repository's cache quota, where
# they would otherwise sit until LRU eviction and push out the caches that are
# still in use (apt, clang-19 .debs, the chr22 fixture). After the save, delete
# every entry under the same prefix on refs/heads/main that was created before
# <saved-key>. Entries created after it (a later push that saved first) are
# left alone, so two overlapping pushes never delete each other's newest entry.
# Without <saved-key> (or if it is not listed, e.g. another job saves it and has
# not yet), the newest entry under the prefix is the one kept.
#
# Needs GH_TOKEN with `actions: write` and GH_REPO. Best effort: a failure to
# list or delete is reported and ignored, since it costs quota, not results.

set -uo pipefail

prefix=${1:?key prefix}
saved=${2:-}

if ! list=$(gh cache list --ref refs/heads/main --key "$prefix" --limit 100 \
    --json id,key,createdAt --jq '.[] | [.id, .key, .createdAt] | @tsv'); then
    echo "::warning::could not list caches under '$prefix'; nothing pruned"
    exit 0
fi

saved_at=$(awk -F'\t' -v k="$saved" '$2 == k { print $3 }' <<< "$list")
if [ -z "$saved_at" ]; then
    # Keep the newest entry instead.
    newest=$(sort -t$'\t' -k3,3 <<< "$list" | tail -n 1)
    saved=$(cut -f2 <<< "$newest")
    saved_at=$(cut -f3 <<< "$newest")
fi
[ -n "$saved_at" ] || {
    echo "no caches under '$prefix'; nothing to prune"
    exit 0
}

pruned=0
while IFS=$'\t' read -r id key created; do
    [ -n "$id" ] || continue
    # ISO-8601 UTC timestamps compare correctly as strings.
    if [ "$key" != "$saved" ] && [[ $created < $saved_at ]]; then
        if gh cache delete "$id" > /dev/null; then
            echo "pruned $key (created $created)"
            pruned=$((pruned + 1))
        else
            echo "::warning::could not delete cache $key"
        fi
    fi
done <<< "$list"
echo "pruned $pruned superseded entr$([ "$pruned" = 1 ] && echo y || echo ies) under '$prefix'"
