#!/bin/bash
###############################################################################
# push_results.sh — Copy results from process/ to git-trackable dir and push
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

set -e
cd /dfs7/adl/tdlong/fly_pool/XQTL2

RESDIR=scripts_oneoffs/ZINC2/debug_chrX_results
mkdir -p ${RESDIR}

# Copy diagnostic file from gitignored process/ to trackable dir
cp -v process/ZINC2/ZINC2_F_v3/ZINC2_F_v3.smooth_diag.chrX.txt \
      ${RESDIR}/smooth_diag.chrX.txt

git pull dev main --rebase
git add ${RESDIR}/
git commit -m 'smooth diag chrX results'
git push dev HEAD:main

echo "Done. Pull from XQTL2-dev to see results."
