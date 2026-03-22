#!/bin/bash
# malathion_test_v2.sh — full scan + figures for the malathion test dataset
#
# Two-replicate experiment: female and male pools treated as biological
# replicates (selected vs control for each sex).  Haplotypes already exist
# in process/malathion_test/, so this starts at Step 5.
#
# Run from XQTL2 project root:
#   bash helpfiles/malathion_test/malathion_test_v2.sh
set -e

DESIGN=helpfiles/malathion_test/design.txt
DIR=process/malathion_test
SCAN=MALATHION_TEST_v2_smooth250
SCAN_DIR=${DIR}/${SCAN}

# ── Haplotype scan (smooth → Wald + H² → concat) ────────────────────────────
scan_out=$(bash scripts/run_scan.sh \
    --design ${DESIGN} \
    --dir    ${DIR} \
    --scan   ${SCAN})
echo "$scan_out"
jid_hap=$(echo "$scan_out" | grep "^done:" | awk '{print $2}')

# ── SNP scan (uses smoothed data from above) ────────────────────────────────
snp_out=$(bash scripts/run_snp_scan.sh \
    --design    ${DESIGN} \
    --dir       ${DIR} \
    --scan      ${SCAN} \
    --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founders  A1,A2,A3,A4,A5,A6,A7,AB8)
echo "$snp_out"
jid_snp=$(echo "$snp_out" | grep "^done:" | awk '{print $2}')

# ── Figures (run after both scans finish) ────────────────────────────────────
sbatch --dependency=afterok:${jid_hap},afterok:${jid_snp} \
    -A tdlong_lab -p standard --mem=8G --time=0:30:00 \
    --wrap="module load R/4.2.2 && \
Rscript scripts/plot_pseudoscan.R \
    --scan ${SCAN_DIR}/${SCAN}.scan.txt \
    --out  ${SCAN_DIR}/${SCAN}.wald.png \
    --format powerpoint --threshold 10 && \
Rscript scripts/plot_H2_overlay.R \
    --scan ${SCAN_DIR}/${SCAN}.scan.txt \
    --out  ${SCAN_DIR}/${SCAN}.H2.png \
    --format powerpoint && \
Rscript scripts/plot_freqsmooth_snp.R \
    --scan ${SCAN_DIR}/${SCAN}.snp_scan.txt \
    --out  ${SCAN_DIR}/${SCAN}.snp.wald.png \
    --format powerpoint --threshold 10"

echo "All jobs submitted."
