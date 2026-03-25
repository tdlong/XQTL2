#!/bin/bash
# Survey repo state: tracked, untracked, and ignored files with sizes
# Run from repo root. Paste output back to Claude.

echo "====== GIT BRANCH / REMOTE STATE ======"
git branch -vv
echo ""
git remote -v
echo ""

echo "====== TRACKED FILES WITH SIZES ======"
git ls-files | xargs -d '\n' du -sh 2>/dev/null | sort -h
echo ""

echo "====== UNTRACKED FILES (not in .gitignore) ======"
git ls-files --others --exclude-standard | xargs -d '\n' du -sh 2>/dev/null | sort -h
echo ""

echo "====== IGNORED FILES (caught by .gitignore) ======"
git ls-files --others --ignored --exclude-standard | xargs -d '\n' du -sh 2>/dev/null | sort -h
echo ""

echo "====== ALL FILES RECURSIVELY WITH SIZES (top 60 largest) ======"
find . -not -path './.git/*' -type f | xargs du -sh 2>/dev/null | sort -rh | head -60
echo ""

echo "====== DIRECTORY SIZES ======"
find . -maxdepth 2 -not -path './.git/*' -type d | xargs du -sh 2>/dev/null | sort -rh
echo ""

echo "====== MODIFIED TRACKED FILES ======"
git diff --name-only
echo ""
