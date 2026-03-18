#!/bin/bash
# fix_move_freqsmooth_stable.sh
# Move freqsmooth and stable scripts from analysis/ back to scripts/
# so they can be tracked in git.

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "Working in: $ROOT"

FILES=(
    haps2scan.freqsmooth.R
    haps2scan.freqsmooth.sh
    haps2scan.stable.R
    haps2scan.stable.sh
    haps2scan.stable.code.R
)

for f in "${FILES[@]}"; do
    src="analysis/$f"
    if [ -f "$src" ]; then
        echo "  mv $src -> scripts/$f"
        mv "$src" "scripts/$f"
    fi
done

echo ""
echo "Done. Run: git add scripts/haps2scan.freqsmooth.* scripts/haps2scan.stable.* && git commit"
