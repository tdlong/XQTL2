#!/bin/bash
# One-time setup: organize public/private files into main and dev branches.
# Run from repo root on the server.
set -e

echo "=== Step 1: Commit public files to main ==="
git checkout main
git add helpfiles/generic_haplotype_parameters.R
git add helpfiles/A_generic_haplotype_parameters.R
git add helpfiles/founder.bams.txt
git add helpfiles/flymap.r6.txt
git add scripts/XQTL_plotting_functions.R
mkdir -p configs && touch configs/.gitkeep
git add configs/.gitkeep
git commit -m "Add helpfiles root files, configs placeholder, XQTL plotting functions"
git push origin main
git push dev main

echo "=== Step 2: Create dev branch ==="
git checkout -b dev

echo "=== Step 3: Rescue scripts buried in ignored directories ==="
# getdata scripts
cp data/raw/AGE_SY/Aug8_25/getdata.sh        scripts_oneoffs/getdata_AGE_SY_Aug8_25.sh
cp data/raw/AGE_SY/May25/getdata_May25.sh    scripts_oneoffs/getdata_AGE_SY_May25.sh
cp data/raw/AGE_SY/Oct28_25/getdata.sh       scripts_oneoffs/getdata_AGE_SY_Oct28_25.sh
cp data/raw/AGE_SY/Sept_2/getdata.sh         scripts_oneoffs/getdata_AGE_SY_Sept_2.sh
cp data/raw/AGE_SY/Sept_24/getdata.sh        scripts_oneoffs/getdata_AGE_SY_Sept_24.sh
cp data/raw/AGE_SY/Sept_24/getdata2.sh       scripts_oneoffs/getdata_AGE_SY_Sept_24b.sh
cp data/raw/JUICE/May25/getdata_May25.sh     scripts_oneoffs/getdata_JUICE_May25.sh
cp data/raw/JUICE/getdata.sh                 scripts_oneoffs/getdata_JUICE.sh
cp data/raw/STARVE/get_data.sh               scripts_oneoffs/getdata_STARVE.sh
cp data/raw/STARVE/get_data.2.sh             scripts_oneoffs/getdata_STARVE_2.sh
cp data/raw/YW/experiment_1/getdata.sh       scripts_oneoffs/getdata_YW_exp1.sh
cp data/raw/ZINC2/get_data.sh                scripts_oneoffs/getdata_ZINC2.sh
cp data/raw/ZINC2/get_data.2.sh              scripts_oneoffs/getdata_ZINC2_2.sh
cp process/ZINC2/ZINC2_F_v3/merge_snp_means.R    scripts_oneoffs/merge_snp_means_ZINC2_F_v3.R
cp process/ZINC2/ZINC2_M_v3/merge_snp_means.R    scripts_oneoffs/merge_snp_means_ZINC2_M_v3.R
cp process/ZINC_Hanson/ZINC_Hanson_v3/merge_snp_means.R scripts_oneoffs/merge_snp_means_ZINC_Hanson_v3.R
cp scripts/temp/mergebams.sh                 scripts_oneoffs/mergebams.sh

echo "=== Step 4: Move root-level one-off scripts ==="
mv finish_v3.sh          scripts_oneoffs/
mv run_ZINC2_F_v3.sh     scripts_oneoffs/
mv run_ZINC2_M_v3.sh     scripts_oneoffs/
mv run_ZINC_Hanson_v3.sh scripts_oneoffs/
mv check_cluster_state.sh scripts_oneoffs/
mv scripts/BU_REFALT2haps.Andreas.code.R scripts_oneoffs/

echo "=== Step 5: Commit all private content to dev ==="
git add scripts_oneoffs/
git add helpfiles/STARVE/ helpfiles/MALATHION/ helpfiles/AGE_SY/
git add helpfiles/ZINC2/ helpfiles/B6885/ helpfiles/MTX/ helpfiles/YW/
git add helpfiles/JUICE/ helpfiles/ZINC_Hanson/ helpfiles/AGE_Aug13_24/
git add helpfiles/pupalHeight2/bam_list.txt
git add configs/
git commit -m "Add private helpfiles, scripts_oneoffs, configs to dev branch"

echo "=== Step 6: Push dev branch to private remote ==="
git push -u dev dev:main

echo ""
echo "Done. Workflow going forward:"
echo "  Public work:  git checkout main  -> git push origin main"
echo "  Private work: git checkout dev   -> git push dev dev:main"
echo "  Keep dev up to date with main:   git checkout dev && git merge main"
