#!/bin/bash
# check_readme.sh — verify that every pipeline/scripts/* and pipeline/helpfiles/*
# path referenced in README.md actually exists on disk.
#
# Usage (from the repo root, or from a project root where pipeline/ is a symlink):
#   bash pipeline/scripts/check_readme.sh
#   bash scripts/check_readme.sh
#
# Exit code: 0 if all references resolve, 1 if any are missing.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
README="$REPO_ROOT/README.md"

if [[ ! -f "$README" ]]; then
    echo "ERROR: README.md not found at $README" >&2
    exit 1
fi

errors=0

check_refs() {
    local pattern="$1"
    while IFS= read -r ref; do
        # pipeline/scripts/foo -> $REPO_ROOT/scripts/foo
        local rel="${ref#pipeline/}"
        local path="$REPO_ROOT/$rel"
        if [[ ! -e "$path" ]]; then
            echo "MISSING: $ref"
            errors=$((errors + 1))
        fi
    done < <(grep -oE "$pattern" "$README" | sort -u)
}

echo "Checking pipeline/scripts/* references..."
check_refs 'pipeline/scripts/[A-Za-z0-9_.]+\.(sh|R)'

echo "Checking pipeline/helpfiles/* references..."
check_refs 'pipeline/helpfiles/[A-Za-z0-9_.]+[A-Za-z0-9_.]'

echo ""
echo "Checking for scripts/* not mentioned in README..."
# Only check .sh files — .R files are implementation details called by their .sh wrappers
# Skip internal/non-user-facing scripts
SKIP="check_readme.sh haps2scan.Apr2025.sh plot_scan_comparison.sh show_project_layout.sh"
while IFS= read -r script; do
    name="$(basename "$script")"
    skip=false
    for s in $SKIP; do [[ "$name" == "$s" ]] && skip=true && break; done
    $skip && continue
    if ! grep -qF "$name" "$README"; then
        echo "UNDOCUMENTED: scripts/$name"
    fi
done < <(find "$REPO_ROOT/scripts" -maxdepth 1 -name "*.sh" | sort)

echo ""
if [[ $errors -gt 0 ]]; then
    echo "FAIL: $errors missing reference(s) found."
    exit 1
else
    echo "OK: all README references resolve."
fi
