#!/bin/bash
# dir_tree.sh — map every directory and its size, save to a file
# Run from anywhere. Pass the root directory as argument, or defaults to current dir.
# Usage: bash dir_tree.sh /dfs7/adl/tdlong/fly_pool > /tmp/dir_tree.txt
#
# Does NOT modify anything.

ROOT=${1:-.}

echo "Directory tree of: $(realpath $ROOT)"
echo "Generated: $(date)"
echo "Host: $(hostname)"
echo ""

# Every directory, sorted, with size
find "$ROOT" -not -path '*/.git/*' -not -name '.git' -type d 2>/dev/null | \
    sort | \
    while read d; do
        size=$(du -sh "$d" 2>/dev/null | cut -f1)
        echo "$size	$d"
    done
