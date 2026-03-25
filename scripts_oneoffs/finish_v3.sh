#!/bin/bash
###############################################################################
# finish_v3.sh — complete the v3 scans
#
# For each scan: merge snp_meansBySample per-chr files, then run figures+tarball.
# The snp_scan.txt already exists; only snp_meansBySample.txt is missing.
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

set -e

finish_scan() {
    local DIR=$1
    local SCAN=$2
    local SCAN_DIR=${DIR}/${SCAN}

    echo "=== ${SCAN} ==="

    # Step 1: merge snp_meansBySample per-chr files using shell cat (no memory issues)
    jid=$(sbatch --parsable \
        -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
        --job-name=snpcat_${SCAN} \
        --wrap="head -1 ${SCAN_DIR}/${SCAN}.snp_meansBySample.chrX.txt > ${SCAN_DIR}/${SCAN}.snp_meansBySample.txt && tail -n +2 -q ${SCAN_DIR}/${SCAN}.snp_meansBySample.chr*.txt >> ${SCAN_DIR}/${SCAN}.snp_meansBySample.txt && echo Done: \$(wc -l < ${SCAN_DIR}/${SCAN}.snp_meansBySample.txt) lines")
    echo "  snp_meansBySample merge: ${jid}"

    # Step 2: figures + tarball (after merge)
    sbatch --dependency=afterok:${jid} \
        -A tdlong_lab -p standard --cpus-per-task=2 --mem-per-cpu=6G --time=1:00:00 \
        --job-name=figs_${SCAN} \
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
    --format powerpoint --threshold 10 && \
cd ${SCAN_DIR} && tar -czf ${SCAN}.tar.gz *.txt *.png"
    echo "  figures + tarball:       submitted (after ${jid})"
}

finish_scan process/ZINC_Hanson ZINC_Hanson_v3
finish_scan process/ZINC2       ZINC2_F_v3
finish_scan process/ZINC2       ZINC2_M_v3

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
