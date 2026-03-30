#!/bin/bash
# server_survey.sh — read-only survey of the full server environment
# Run from the XQTL2 repo root. Paste output back to Claude.
# Does NOT modify anything.

REPO_ROOT=$(pwd)
PARENT=$(dirname "$REPO_ROOT")

echo "======================================================================"
echo "  SERVER SURVEY — $(date)"
echo "  Repo root: $REPO_ROOT"
echo "  Parent:    $PARENT"
echo "  Host:      $(hostname)"
echo "======================================================================"
echo ""

echo "====== GIT: BRANCH / REMOTES ======"
git branch -vv
git remote -v
echo ""

echo "====== GIT: RECENT COMMITS (20) ======"
git log --oneline -20
echo ""

echo "====== GIT: TRACKED FILES MODIFIED (uncommitted) ======"
git diff --name-only
if [ -z "$(git diff --name-only)" ]; then echo "(none)"; fi
echo ""

echo "====== GIT: STAGED BUT NOT COMMITTED ======"
git diff --cached --name-only
if [ -z "$(git diff --cached --name-only)" ]; then echo "(none)"; fi
echo ""

echo "====== GIT: UNTRACKED (not gitignored) ======"
git ls-files --others --exclude-standard
echo ""

echo "====== PARENT DIRECTORY CONTENTS ======"
ls -la "$PARENT/"
echo ""

echo "====== PARENT DIRECTORY: TOP-LEVEL SIZES ======"
du -sh "$PARENT"/*/  2>/dev/null | sort -rh
echo ""

echo "====== REPO ROOT: DIRECTORY SIZES (depth 2) ======"
find . -maxdepth 2 -not -path './.git/*' -type d | \
    xargs du -sh 2>/dev/null | sort -rh
echo ""

echo "====== data/ CONTENTS ======"
if [ -d data ]; then
    find data -maxdepth 3 -type d | head -60
    echo "--- sizes ---"
    du -sh data/*/  2>/dev/null | sort -rh
else
    echo "(no data/ directory)"
fi
echo ""

echo "====== process/ CONTENTS ======"
if [ -d process ]; then
    find process -maxdepth 3 -type d | head -60
    echo "--- sizes ---"
    du -sh process/*/  2>/dev/null | sort -rh
else
    echo "(no process/ directory)"
fi
echo ""

echo "====== scripts_oneoffs/ CONTENTS ======"
if [ -d scripts_oneoffs ]; then
    find scripts_oneoffs -maxdepth 3 -type f | sort
else
    echo "(no scripts_oneoffs/ directory)"
fi
echo ""

echo "====== helpfiles/ CONTENTS ======"
if [ -d helpfiles ]; then
    find helpfiles -maxdepth 3 -type f | sort
else
    echo "(no helpfiles/ directory)"
fi
echo ""

echo "====== figures/ CONTENTS ======"
if [ -d figures ]; then
    find figures -maxdepth 3 -type f | sort | head -60
    echo "--- sizes ---"
    du -sh figures/ 2>/dev/null
else
    echo "(no figures/ directory)"
fi
echo ""

echo "====== output/ CONTENTS ======"
if [ -d output ]; then
    find output -maxdepth 3 -type f | sort | head -60
    du -sh output/ 2>/dev/null
else
    echo "(no output/ directory)"
fi
echo ""

echo "====== TOP 30 LARGEST FILES (excluding .git) ======"
find . -not -path './.git/*' -type f -exec du -sh {} + 2>/dev/null | \
    sort -rh | head -30
echo ""

echo "====== .gitignore CONTENTS ======"
if [ -f .gitignore ]; then
    cat .gitignore
else
    echo "(no .gitignore)"
fi
echo ""

echo "====== SURVEY COMPLETE ======"
