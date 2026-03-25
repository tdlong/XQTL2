#!/bin/bash
###############################################################################
# run_ZINC_Hanson_v3.sh — Re-analyze Hanson (pA) data with updated XQTL2 pipeline
#
# Starts from existing haplotype files (R.haps.*.out.rds in process/ZINC_Hanson/).
# Runs haplotype scan (Wald + H2) and SNP scan with 250 kb smoothing (default).
#
# INPUT:  process/ZINC_Hanson/R.haps.*.out.rds  (existing)
#         helpfiles/ZINC_Hanson/ZINC_Hanson.test.txt
#         helpfiles/FREQ_SNPs_Apop.cM.txt.gz
# OUTPUT: process/ZINC_Hanson/ZINC_Hanson_v3/
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
#
# DEPLOYMENT:
#   scp pipeline/run_ZINC_Hanson_v3.sh tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/
#   ssh hpc3 "cd /dfs7/adl/tdlong/fly_pool/XQTL2 && bash run_ZINC_Hanson_v3.sh"
###############################################################################

set -e

DESIGN=helpfiles/ZINC_Hanson/ZINC_Hanson.test.txt
DIR=process/ZINC_Hanson
SCAN=ZINC_Hanson_v3
SNP_TABLE=helpfiles/FREQ_SNPs_Apop.cM.txt.gz
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8

echo "=== ZINC_Hanson reanalysis with updated XQTL2 pipeline ==="
echo "Design:    $DESIGN"
echo "Dir:       $DIR"
echo "Scan:      $SCAN"
echo "SNP table: $SNP_TABLE"
echo ""

# ── Step 5a: haplotype scan (smooth → Wald + H2 → concat) ──────────────────
scan_out=$(bash scripts/run_scan.sh \
    --design       ${DESIGN} \
    --dir          ${DIR} \
    --scan         ${SCAN} \
    --cpus-per-task 2 \
    --mem-per-cpu  6G)
echo "$scan_out"
jid_hap=$(echo "$scan_out" | grep "^done:" | awk '{print $2}')

# ── Step 5b: SNP scan ───────────────────────────────────────────────────────
snp_out=$(bash scripts/run_snp_scan.sh \
    --design    ${DESIGN} \
    --dir       ${DIR} \
    --scan      ${SCAN} \
    --snp-table ${SNP_TABLE} \
    --founders  ${FOUNDERS} \
    --after     ${jid_hap})
echo "$snp_out"
jid_snp=$(echo "$snp_out" | grep "^done:" | awk '{print $2}')

# ── Step 6: figures + tarball ───────────────────────────────────────────────
SCAN_DIR=${DIR}/${SCAN}
sbatch --dependency=afterok:${jid_hap},afterok:${jid_snp} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
    --job-name=figs_ZINC_Hanson_v3 \
    --wrap="module load R/4.2.2 && \
Rscript scripts/plot_pseudoscan.R \
    --scan      ${SCAN_DIR}/${SCAN}.scan.txt \
    --out       ${SCAN_DIR}/${SCAN}.wald.png \
    --format    powerpoint \
    --threshold 10 && \
Rscript scripts/plot_H2_overlay.R \
    --scan   ${SCAN_DIR}/${SCAN}.scan.txt \
    --out    ${SCAN_DIR}/${SCAN}.H2.png \
    --format powerpoint && \
Rscript scripts/plot_freqsmooth_snp.R \
    --scan      ${SCAN_DIR}/${SCAN}.snp_scan.txt \
    --out       ${SCAN_DIR}/${SCAN}.snp.wald.png \
    --format    powerpoint \
    --threshold 10 && \
cd ${SCAN_DIR} && tar -czf ${SCAN}.tar.gz *.txt *.png"

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo "Results: ${SCAN_DIR}/"
echo "Download when done: scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/${SCAN_DIR}/${SCAN}.tar.gz notes/"
