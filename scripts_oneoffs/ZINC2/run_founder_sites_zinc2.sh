#!/bin/bash
###############################################################################
# run_founder_sites_zinc2.sh — Validate founder-sites RefAlt approach on ZINC2
#
# Builds a B-founder site catalog (if not already done), re-derives RefAlt
# for all ZINC2 pools using that catalog, then compares against the original
# joint-called RefAlt in process/ZINC2/.
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
#   bash scripts_oneoffs/ZINC2/run_founder_sites_zinc2.sh
###############################################################################

set -e

BAMS=helpfiles/ZINC2/ZINC2.bams
FOUNDER_SITES=process/B_founder_sites
FROMSITES=process/ZINC2_fromsites
ORIGINAL=process/ZINC2
LOGS=logs/ZINC2_fromsites

mkdir -p $FOUNDER_SITES $FROMSITES $LOGS

echo "=== Founder-sites RefAlt validation — ZINC2 ==="
echo ""

# ── Phase 1: Build founder site catalog (skip if already done) ───────────────
if ls ${FOUNDER_SITES}/founder_sites.chr*.vcf.gz 1>/dev/null 2>&1; then
  echo "Phase 1: founder site catalog already exists in ${FOUNDER_SITES}/ — skipping"
  jid1=""
  dep1=""
else
  echo "Phase 1: Building B-founder site catalog..."
  jid1=$(sbatch scripts_oneoffs/bam2founder_sites.sh \
           helpfiles/B_founders.bams.txt $FOUNDER_SITES \
         | awk '{print $4}')
  echo "  Submitted job $jid1 (array 1-5)"
  echo -e "founderSites\t$jid1" >> $LOGS/jobs.log
  dep1="--dependency=afterok:$jid1"
fi

# ── Phase 2: Per-sample RefAlt at founder sites ───────────────────────────────
echo "Phase 2: Submitting per-sample jobs..."
jid2_dep=""
n=0
while read bam; do
  [[ $bam == /dfs7* ]] && continue   # skip founder lines
  sample=$(basename $bam .bam)
  jid=$(sbatch $dep1 \
               --output=${LOGS}/${sample}.%A_%a.out \
               scripts_oneoffs/bam2bcf2REFALT_fromsites.sh \
               $bam $FOUNDER_SITES $FROMSITES \
         | awk '{print $4}')
  jid2_dep="${jid2_dep:+$jid2_dep:}$jid"
  echo -e "${sample}\t$jid" >> $LOGS/jobs.log
  n=$((n+1))
done < $BAMS
echo "  Submitted $n sample jobs"

# ── Phase 3: Consolidate per-sample files ────────────────────────────────────
echo "Phase 3: Consolidating per-sample RefAlt files..."
jid3=$(sbatch --dependency=afterok:$jid2_dep \
              --output=${LOGS}/consolidate.%A_%a.out \
              scripts_oneoffs/consolidate_refalt_fromsites.sh \
              $FROMSITES $FOUNDER_SITES \
        | awk '{print $4}')
echo "  Submitted job $jid3 (array 1-5)"
echo -e "consolidate\t$jid3" >> $LOGS/jobs.log

# ── Phase 4: Compare against original ZINC2 RefAlt ───────────────────────────
echo "Phase 4: Comparing against original RefAlt..."
jid4=$(sbatch --dependency=afterok:$jid3 \
              --output=${LOGS}/compare.%A.out \
              scripts_oneoffs/compare_refalt_zinc2.sh \
              $ORIGINAL $FROMSITES \
        | awk '{print $4}')
echo "  Submitted job $jid4"
echo -e "compare\t$jid4" >> $LOGS/jobs.log

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo "Comparison results: cat ${LOGS}/compare.${jid4}.out"
echo ""
echo "After jobs finish, profile resource usage:"
echo "  bash scripts_oneoffs/profile_jobs.sh ${LOGS}/jobs.log"
