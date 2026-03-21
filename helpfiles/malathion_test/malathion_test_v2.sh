#!/bin/bash
# malathion_test_v2.sh
#
# Test the freqsmooth pipeline on the malathion dataset.
# Run from the XQTL2 project root.
#
# Prerequisites:
#   1. R.haps.<chr>.out.rds files exist in process/malathion_test/
#      (steps 1-4 of malathion_test_pipeline.sh already completed)
#   2. SNP_TABLE path below points to FREQ_SNPs.cM.txt.gz on this cluster
#   3. helpfiles/malathion_test/MALATHION_TEST_v2_smooth125.R figure script exists
#
# Submit with:  bash helpfiles/malathion_test/malathion_test_v2.sh

PROJECT=malathion_test
DESIGN=helpfiles/${PROJECT}/design.txt
SCAN=MALATHION_TEST_v2_smooth125
SMOOTH_KB=125
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8
FIGURE=helpfiles/${PROJECT}/MALATHION_TEST_v2_smooth125.R

# ── SET THIS to the cluster path of FREQ_SNPs.cM.txt.gz ──────────────────────
SNP_TABLE=/dfs7/adl/tdlong/fly_pool/FREQ_SNPs_Apop.cM.txt.gz

# ── Step 5a: smooth haplotype frequencies ─────────────────────────────────────
# Reads R.haps.<chr>.out.rds, writes smoothed RDS + meansBySample per chromosome
jid_smooth=$(sbatch --parsable \
    --array=1-5 scripts_freqsmooth/smooth_haps.sh \
    --rfile     ${DESIGN}            \
    --dir       process/${PROJECT}   \
    --outdir    ${SCAN}              \
    --smooth-kb ${SMOOTH_KB})
echo "smooth:  $jid_smooth"

# Explicit per-task dependency: afterok on an array job ID is unreliable in
# older SLURM (may release on first task completion).  List all 5 task IDs.
smooth_dep="afterok:${jid_smooth}_1:${jid_smooth}_2:${jid_smooth}_3:${jid_smooth}_4:${jid_smooth}_5"

# ── Step 5b: Wald scan on smoothed frequencies ────────────────────────────────
# Reads smoothed RDS, writes scan.<chr>.txt per chromosome
jid_scan=$(sbatch --parsable --dependency=${smooth_dep} \
    --array=1-5 scripts_freqsmooth/freqsmooth_scan.sh \
    --rfile   ${DESIGN}            \
    --dir     process/${PROJECT}   \
    --outdir  ${SCAN})
echo "scan:    $jid_scan"

# ── Step 5c: SNP scan (runs in parallel with 5b) ──────────────────────────────
# Reads smoothed RDS + FREQ_SNPs, writes snp_scan.<chr>.txt per chromosome
jid_snp=$(sbatch --parsable --dependency=${smooth_dep} \
    --array=1-5 scripts_freqsmooth/snp_scan.sh \
    --rfile     ${DESIGN}            \
    --dir       process/${PROJECT}   \
    --outdir    ${SCAN}              \
    --snp-table ${SNP_TABLE}         \
    --founders  ${FOUNDERS})
echo "snp:     $jid_snp"

# ── Step 6: concatenate chromosomes ──────────────────────────────────────────
scan_dep="afterok:${jid_scan}_1:${jid_scan}_2:${jid_scan}_3:${jid_scan}_4:${jid_scan}_5"
snp_dep="afterok:${jid_snp}_1:${jid_snp}_2:${jid_snp}_3:${jid_snp}_4:${jid_snp}_5"
jid_concat=$(sbatch --parsable --dependency=${scan_dep}:${snp_dep#afterok:} \
    -A tdlong_lab -p standard --mem=10G \
    --wrap="bash scripts/concat_Chromosome_Scans.sh process/${PROJECT}/${SCAN}")
echo "concat:  $jid_concat"

# ── Step 7: publication figure ────────────────────────────────────────────────
jid_fig=$(sbatch --parsable --dependency=afterok:${jid_concat} \
    -A tdlong_lab -p standard --mem=10G \
    --wrap="module load R/4.2.2 && Rscript ${FIGURE}")
echo "figure:  $jid_fig"
