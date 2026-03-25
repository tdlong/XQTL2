#!/bin/bash
# Summarize repo state. Run from repo root. Paste output back to Claude.

echo "====== BRANCH / REMOTES ======"
git branch -vv
git remote -v
echo ""

echo "====== MODIFIED TRACKED FILES ======"
git diff --name-only
echo ""

echo "====== UNTRACKED (not gitignored) — directory sizes ======"
git ls-files --others --exclude-standard | \
    awk -F/ '{print $1}' | sort -u | \
    while read d; do du -sh "$d" 2>/dev/null; done
echo ""

echo "====== IGNORED — top 30 largest files ======"
git ls-files --others --ignored --exclude-standard | \
    xargs du -sh 2>/dev/null | sort -rh | head -30
echo ""

echo "====== IGNORED — file extensions summary ======"
git ls-files --others --ignored --exclude-standard | \
    grep -o '\.[^./]*$' | sort | uniq -c | sort -rn | head -20
echo ""

echo "====== TOP 20 LARGEST FILES (anywhere, not .git) ======"
find . -not -path './.git/*' -type f -exec du -sh {} + 2>/dev/null | \
    sort -rh | head -20
echo ""

echo "====== DIRECTORY SIZES (depth 2) ======"
find . -maxdepth 2 -not -path './.git/*' -type d | \
    xargs du -sh 2>/dev/null | sort -rh
echo ""

echo "====== IGNORED .sh AND .R FILES (potential lost scripts) ======"
git ls-files --others --ignored --exclude-standard | grep -E '\.(sh|R)$'
echo ""
