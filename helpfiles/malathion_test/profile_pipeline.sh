#!/bin/bash
# profile_pipeline.sh — run malathion scan pipeline and collect seff profiling
#
# Runs the full scan (smooth → hap_scan → concat → snp_scan → snp_concat →
# figures), then submits a final job that runs seff on every job and writes
# results to a single file.
#
# Run from XQTL2 project root:
#   bash helpfiles/malathion_test/profile_pipeline.sh
#
# Output: process/malathion_test/MALATHION_TEST_v2_smooth250/seff_profile.txt
set -e

DESIGN=helpfiles/malathion_test/design.txt
DIR=process/malathion_test
SCAN=MALATHION_TEST_v2_smooth250
SCAN_DIR=${DIR}/${SCAN}
OUTFILE=${SCAN_DIR}/seff_profile.txt

# ── Haplotype scan ───────────────────────────────────────────────────────────
scan_out=$(bash scripts/run_scan.sh \
    --design ${DESIGN} \
    --dir    ${DIR} \
    --scan   ${SCAN})
echo "$scan_out"

jid_smooth=$(echo "$scan_out" | grep "^smooth:" | awk '{print $2}')
jid_hapscan=$(echo "$scan_out" | grep "^hap_scan:" | awk '{print $2}')
jid_concat=$(echo "$scan_out" | grep "^concat:" | awk '{print $2}')

# ── SNP scan ─────────────────────────────────────────────────────────────────
snp_out=$(bash scripts/run_snp_scan.sh \
    --design    ${DESIGN} \
    --dir       ${DIR} \
    --scan      ${SCAN} \
    --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founders  A1,A2,A3,A4,A5,A6,A7,AB8)
echo "$snp_out"

jid_snpscan=$(echo "$snp_out" | grep "^snp_scan:" | awk '{print $2}')
jid_snpconcat=$(echo "$snp_out" | grep "^snp_concat:" | awk '{print $2}')

# ── Figures ──────────────────────────────────────────────────────────────────
jid_figures=$(sbatch --parsable \
    --dependency=afterok:${jid_concat},afterok:${jid_snpconcat} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=0:10:00 \
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
cd ${SCAN_DIR} && tar -czf ${SCAN}.tar.gz *.txt *.png")
echo "figures:    $jid_figures"

# ── Collect seff profiling after everything finishes ─────────────────────────
ALL_JIDS="${jid_smooth},${jid_hapscan},${jid_concat},${jid_snpscan},${jid_snpconcat},${jid_figures}"

jid_seff=$(sbatch --parsable \
    --dependency=afterany:${jid_figures} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=1G --time=0:10:00 \
    --wrap="
echo '=== SLURM Resource Profiling ===' > ${OUTFILE}
echo 'Date:' \$(date) >> ${OUTFILE}
echo 'Pipeline: malathion_test_v2 (from haps)' >> ${OUTFILE}
echo '' >> ${OUTFILE}

echo '--- smooth_haps (array job ${jid_smooth}) ---' >> ${OUTFILE}
seff ${jid_smooth}_1 >> ${OUTFILE} 2>&1
echo '' >> ${OUTFILE}

echo '--- hap_scan (array job ${jid_hapscan}) ---' >> ${OUTFILE}
seff ${jid_hapscan}_1 >> ${OUTFILE} 2>&1
echo '' >> ${OUTFILE}

echo '--- concat_scans (job ${jid_concat}) ---' >> ${OUTFILE}
seff ${jid_concat} >> ${OUTFILE} 2>&1
echo '' >> ${OUTFILE}

echo '--- snp_scan (array job ${jid_snpscan}) ---' >> ${OUTFILE}
seff ${jid_snpscan}_1 >> ${OUTFILE} 2>&1
echo '' >> ${OUTFILE}

echo '--- snp_concat (job ${jid_snpconcat}) ---' >> ${OUTFILE}
seff ${jid_snpconcat} >> ${OUTFILE} 2>&1
echo '' >> ${OUTFILE}

echo '--- figures (job ${jid_figures}) ---' >> ${OUTFILE}
seff ${jid_figures} >> ${OUTFILE} 2>&1
echo '' >> ${OUTFILE}

echo '=== All job IDs ===' >> ${OUTFILE}
echo 'smooth:     ${jid_smooth}' >> ${OUTFILE}
echo 'hap_scan:   ${jid_hapscan}' >> ${OUTFILE}
echo 'concat:     ${jid_concat}' >> ${OUTFILE}
echo 'snp_scan:   ${jid_snpscan}' >> ${OUTFILE}
echo 'snp_concat: ${jid_snpconcat}' >> ${OUTFILE}
echo 'figures:    ${jid_figures}' >> ${OUTFILE}
echo 'Done.' >> ${OUTFILE}
")
echo "seff:       $jid_seff"

echo ""
echo "All jobs submitted. Job IDs:"
echo "  smooth:     ${jid_smooth}"
echo "  hap_scan:   ${jid_hapscan}"
echo "  concat:     ${jid_concat}"
echo "  snp_scan:   ${jid_snpscan}"
echo "  snp_concat: ${jid_snpconcat}"
echo "  figures:    ${jid_figures}"
echo "  seff:       ${jid_seff}"
echo ""
echo "When done, results in: ${OUTFILE}"
