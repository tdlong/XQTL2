#!/bin/bash
# malathion_test_v2.sh
#
# Test the freqsmooth pipeline on the malathion dataset.
# Run from the XQTL2 project root.
#
# Prerequisites:
#   1. R.haps.<chr>.out.rds files exist in process/malathion_test/
#      (steps 1-4 of malathion_test_pipeline.sh already completed)
#   2. helpfiles/malathion_test/MALATHION_TEST_v2_smooth125.R figure script exists
#
# Submit with:  bash helpfiles/malathion_test/malathion_test_v2.sh

set -e   # exit immediately if any sbatch fails

PROJECT=malathion_test
DESIGN=helpfiles/${PROJECT}/design.txt
SCAN=MALATHION_TEST_v2_smooth125
SMOOTH_KB=125
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8
SNP_TABLE=helpfiles/FREQ_SNPs_Apop.cM.txt.gz
FIGURE=helpfiles/${PROJECT}/MALATHION_TEST_v2_smooth125.R

# ── Step 5a: smooth haplotype frequencies ─────────────────────────────────────
jid_smooth=$(sbatch --parsable \
    --array=1-5 scripts/smooth_haps.sh \
    --rfile     ${DESIGN}          \
    --dir       process/${PROJECT} \
    --outdir    ${SCAN}            \
    --smooth-kb ${SMOOTH_KB})
echo "smooth:  $jid_smooth"

# ── Step 5b: haplotype Wald scan ──────────────────────────────────────────────
jid_scan=$(sbatch --parsable --dependency=afterok:${jid_smooth} \
    --array=1-5 scripts/hap_scan.sh \
    --rfile   ${DESIGN}          \
    --dir     process/${PROJECT} \
    --outdir  ${SCAN})
echo "scan:    $jid_scan"

# ── Step 5c: SNP scan (add-on; runs in parallel with 5b) ──────────────────────
jid_snp=$(sbatch --parsable --dependency=afterok:${jid_smooth} \
    --array=1-5 scripts/snp_scan.sh \
    --rfile     ${DESIGN}          \
    --dir       process/${PROJECT} \
    --outdir    ${SCAN}            \
    --snp-table ${SNP_TABLE}       \
    --founders  ${FOUNDERS})
echo "snp:     $jid_snp"

# ── Step 6: concatenate haplotype scan chromosomes ────────────────────────────
# Depends only on the haplotype scan — SNP scan is independent
jid_concat=$(sbatch --parsable --dependency=afterok:${jid_scan} \
    -A tdlong_lab -p standard \
    --cpus-per-task=2 --mem-per-cpu=6G \
    --wrap="bash scripts/concat_Chromosome_Scans.sh process/${PROJECT}/${SCAN}")
echo "concat:  $jid_concat"

# ── Step 7: publication figure (Wald + Falconer H2 + Cutler H2) ───────────────
jid_fig=$(sbatch --parsable --dependency=afterok:${jid_concat} \
    -A tdlong_lab -p standard \
    --cpus-per-task=2 --mem-per-cpu=6G \
    --wrap="module load R/4.2.2 && Rscript ${FIGURE}")
echo "figure:  $jid_fig"
