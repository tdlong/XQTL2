#!/bin/bash
###############################################################################
# resume_v3.sh — Resume v3 scans after concat_scans.sh --snp OOM fix (issue #3)
#
# Context: rerun_v3.sh completed smooth_haps, hap_scan, hap_concat, and
# snp_scan for all three scans. snp_concat OOMed (issue #3). This script
# picks up from there: snp_concat → figs + tarball → snp_comparison.
#
# INPUT:  process/ZINC_Hanson/ZINC_Hanson_v3/snp_scan.chr*.txt  (existing)
#         process/ZINC2/ZINC2_F_v3/snp_scan.chr*.txt            (existing)
#         process/ZINC2/ZINC2_M_v3/snp_scan.chr*.txt            (existing)
# OUTPUT: process/ZINC_Hanson/ZINC_Hanson_v3/ZINC_Hanson_v3.tar.gz
#         process/ZINC2/ZINC2_F_v3/ZINC2_F_v3.tar.gz
#         process/ZINC2/ZINC2_M_v3/ZINC2_M_v3.tar.gz
#         output/snp_shared_poly.tsv
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

set -e

WD=/dfs7/adl/tdlong/fly_pool/XQTL2

# ── Step 0: pull latest pipeline (issue #3 fix) ───────────────────────────────
echo "=== Pulling latest XQTL2 pipeline ==="
git pull
git -C scripts pull
echo ""

# ── Helper: snp_concat + figs for one scan, returns figs job ID ──────────────
resume_scan() {
    local DIR=$1
    local SCAN=$2

    local SCAN_DIR=${DIR}/${SCAN}
    echo "=== ${SCAN} ==="

    local jid_concat=$(sbatch --parsable \
        -A tdlong_lab -p standard --cpus-per-task=2 --mem-per-cpu=6G --time=1:00:00 \
        --job-name=snp_concat_${SCAN} \
        --wrap="bash scripts/concat_scans.sh --snp ${SCAN_DIR}")
    echo "  snp_concat: ${jid_concat}"

    local jid_figs=$(sbatch --parsable \
        --dependency=afterok:${jid_concat} \
        -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
        --job-name=figs_${SCAN} \
        --wrap="module load R/4.2.2 && \
Rscript scripts/plot_pseudoscan.R --scan ${SCAN_DIR}/${SCAN}.scan.txt --out ${SCAN_DIR}/${SCAN}.wald.png --format powerpoint --threshold 10 && \
Rscript scripts/plot_H2_overlay.R --scan ${SCAN_DIR}/${SCAN}.scan.txt --out ${SCAN_DIR}/${SCAN}.H2.png --format powerpoint && \
Rscript scripts/plot_freqsmooth_snp.R --scan ${SCAN_DIR}/${SCAN}.snp_scan.txt --out ${SCAN_DIR}/${SCAN}.snp.wald.png --format powerpoint --threshold 10 && \
cd ${SCAN_DIR} && tar -czf ${SCAN}.tar.gz *.txt *.png")
    echo "  figs + tarball: ${jid_figs}"

    echo "done_figs: ${jid_figs}"
}

# ── Submit snp_concat + figs for all three scans ─────────────────────────────
out_hanson=$(resume_scan process/ZINC_Hanson ZINC_Hanson_v3)
echo "$out_hanson"
jid_figs_hanson=$(echo "$out_hanson" | grep "^done_figs:" | awk '{print $2}')

out_F=$(resume_scan process/ZINC2 ZINC2_F_v3)
echo "$out_F"
jid_figs_F=$(echo "$out_F" | grep "^done_figs:" | awk '{print $2}')

out_M=$(resume_scan process/ZINC2 ZINC2_M_v3)
echo "$out_M"

# ── SNP comparison: after Hanson and ZINC2_F figs complete ───────────────────
echo ""
echo "=== SNP comparison (server-side) ==="
mkdir -p ${WD}/output
jid_snpcomp=$(sbatch --parsable \
    --dependency=afterok:${jid_figs_hanson},afterok:${jid_figs_F} \
    -A tdlong_lab -p standard --cpus-per-task=2 --mem-per-cpu=6G --time=2:00:00 \
    --job-name=snp_comparison \
    --wrap="module load R/4.2.2 && cd ${WD} && Rscript scripts_oneoffs/ZINC2/snp_comparison_server.R")
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
