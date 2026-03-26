#!/bin/bash
###############################################################################
# rerun_v3.sh — Rerun v3 scans after XQTL2 bug fix (issue #2; overwrites v3)
#
# Bug fixed: per-founder masking before smoother; unresolvable windows no
# longer produce out-of-[0,1] frequencies (chrX ~21.5 Mb spike).
#
# Steps:
#   1. git pull the XQTL2 pipeline scripts
#   2. Submit ZINC_Hanson_v3, ZINC2_F_v3, ZINC2_M_v3 scans (overwrites existing)
#      (each: smooth → hap scan → SNP scan → figs + tarball)
#   3. After ZINC_Hanson_v3 and ZINC2_F_v3 figs complete:
#      run snp_comparison_server.R to produce output/snp_shared_poly.tsv
#
# INPUT:  process/ZINC_Hanson/R.haps.*.out.rds  (existing)
#         process/ZINC2/R.haps.*.out.rds          (existing)
# OUTPUT: process/ZINC_Hanson/ZINC_Hanson_v3/    (overwritten)
#         process/ZINC2/ZINC2_F_v3/              (overwritten)
#         process/ZINC2/ZINC2_M_v3/              (overwritten)
#         output/snp_shared_poly.tsv
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
#
# DEPLOYMENT:
#   Push to XQTL2-dev, then on cluster:
#   ssh hpc3 "cd /dfs7/adl/tdlong/fly_pool/XQTL2 && bash scripts_oneoffs/Zinc2/rerun_v3.sh"
###############################################################################

set -e

WD=/dfs7/adl/tdlong/fly_pool/XQTL2

# ── Step 0: pull latest pipeline and dev scripts ─────────────────────────────
echo "=== Pulling latest XQTL2 pipeline ==="
git pull
git -C scripts pull
echo ""

# ── Helper: submit one scan, return figs job ID ──────────────────────────────
submit_scan() {
    local DESIGN=$1
    local DIR=$2
    local SCAN=$3
    local SNP_TABLE=$4
    local FOUNDERS=$5

    echo "=== ${SCAN} ==="

    scan_out=$(bash scripts/run_scan.sh \
        --design        ${DESIGN} \
        --dir           ${DIR} \
        --scan          ${SCAN} \
        --cpus-per-task 2 \
        --mem-per-cpu   6G)
    echo "$scan_out"
    local jid_hap=$(echo "$scan_out" | grep "^done:" | awk '{print $2}')

    snp_out=$(bash scripts/run_snp_scan.sh \
        --design    ${DESIGN} \
        --dir       ${DIR} \
        --scan      ${SCAN} \
        --snp-table ${SNP_TABLE} \
        --founders  ${FOUNDERS} \
        --after     ${jid_hap})
    echo "$snp_out"
    local jid_snp=$(echo "$snp_out" | grep "^done:" | awk '{print $2}')

    local SCAN_DIR=${DIR}/${SCAN}
    local jid_figs=$(sbatch --parsable \
        --dependency=afterok:${jid_hap},afterok:${jid_snp} \
        -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
        --job-name=figs_${SCAN} \
        --wrap="module load R/4.2.2 && Rscript scripts/plot_pseudoscan.R --scan ${SCAN_DIR}/${SCAN}.scan.txt --out ${SCAN_DIR}/${SCAN}.wald.png --format powerpoint --threshold 10 && Rscript scripts/plot_H2_overlay.R --scan ${SCAN_DIR}/${SCAN}.scan.txt --out ${SCAN_DIR}/${SCAN}.H2.png --format powerpoint && Rscript scripts/plot_freqsmooth_snp.R --scan ${SCAN_DIR}/${SCAN}.snp_scan.txt --out ${SCAN_DIR}/${SCAN}.snp.wald.png --format powerpoint --threshold 10 && cd ${SCAN_DIR} && tar -czf ${SCAN}.tar.gz *.txt *.png")
    echo "  figs + tarball: ${jid_figs}"

    echo "done_figs: ${jid_figs}"
}

# ── Submit all three scans ───────────────────────────────────────────────────
out_hanson=$(submit_scan \
    helpfiles/ZINC_Hanson/ZINC_Hanson.test.txt \
    process/ZINC_Hanson \
    ZINC_Hanson_v3 \
    helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    A1,A2,A3,A4,A5,A6,A7,AB8)
jid_figs_hanson=$(echo "$out_hanson" | grep "^done_figs:" | awk '{print $2}')

out_F=$(submit_scan \
    helpfiles/ZINC2/Zinc2.test.F.txt \
    process/ZINC2 \
    ZINC2_F_v3 \
    helpfiles/FREQ_SNPs_Bpop.cM.txt.gz \
    B1,B2,B3,B4,B5,B6,B7,AB8)
jid_figs_F=$(echo "$out_F" | grep "^done_figs:" | awk '{print $2}')

submit_scan \
    helpfiles/ZINC2/Zinc2.test.M.txt \
    process/ZINC2 \
    ZINC2_M_v3 \
    helpfiles/FREQ_SNPs_Bpop.cM.txt.gz \
    B1,B2,B3,B4,B5,B6,B7,AB8

# ── SNP comparison: server-side, after Hanson and ZINC2_F figs complete ──────
echo ""
echo "=== SNP comparison (server-side) ==="
mkdir -p ${WD}/output
jid_snpcomp=$(sbatch --parsable \
    --dependency=afterok:${jid_figs_hanson},afterok:${jid_figs_F} \
    -A tdlong_lab -p standard --cpus-per-task=2 --mem-per-cpu=6G --time=2:00:00 \
    --job-name=snp_comparison \
    --wrap="module load R/4.2.2 && cd ${WD} && Rscript scripts_oneoffs/Zinc2/snp_comparison_server.R")
echo "  snp_comparison: ${jid_snpcomp}"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo ""
echo "When complete, download:"
echo "  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC_Hanson/ZINC_Hanson_v3/ZINC_Hanson_v3.tar.gz notes/"
echo "  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/ZINC2_F_v3/ZINC2_F_v3.tar.gz notes/"
echo "  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/process/ZINC2/ZINC2_M_v3/ZINC2_M_v3.tar.gz notes/"
echo "  scp tdlong@hpc3.rcic.uci.edu:/dfs7/adl/tdlong/fly_pool/XQTL2/output/snp_shared_poly.tsv notes/"
