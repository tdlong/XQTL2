  #!/bin/bash
  # malathion_test_pipeline.sh
  # Full pipeline for malathion_test, starting from bam files.
  # Run from the XQTL2 project root.
  #
  # Before running this script:
  #   mkdir -p data/bam/malathion_test
  #   cd data/bam/malathion_test
  #   ln -s /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/data/bam/c_F_1.bam .
  #   ln -s /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/data/bam/c_F_1.bam.bai .
  #   ln -s /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/data/bam/c_M_1.bam .
  #   ln -s /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/data/bam/c_M_1.bam.bai .
  #   ln -s /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/data/bam/s_F_1.bam .
  #   ln -s /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/data/bam/s_F_1.bam.bai .
  #   ln -s /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/data/bam/s_M_1.bam .
  #   ln -s /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/data/bam/s_M_1.bam.bai .
  #   cd ../../..

  PROJECT=malathion_test
  PARFILE=helpfiles/${PROJECT}/hap_params.R
  DESIGN=helpfiles/${PROJECT}/design.txt
  SCAN=MALATHION_TEST_smooth125
  COV_KB=125
  FREQ_KB=125
  FIGURE=helpfiles/${PROJECT}/MALATHION_TEST_smooth125.R

  # ── Build bams file ───────────────────────────────────────────────────────────
  mkdir -p process/${PROJECT}
  find -L data/bam/${PROJECT} -name "*.bam" > helpfiles/${PROJECT}/bams
  cat helpfiles/founder.bams.txt >> helpfiles/${PROJECT}/bams

  # ── Step 3: REFALT counts ─────────────────────────────────────────────────────
  jid_refalt=$(sbatch --parsable \
      --array=1-5 scripts/bam2bcf2REFALT.sh \
      helpfiles/${PROJECT}/bams process/${PROJECT})
  echo "REFALT:  $jid_refalt"

  # ── Step 4: haplotypes ────────────────────────────────────────────────────────
  jid_haps=$(sbatch --parsable --dependency=afterok:${jid_refalt} \
      --array=1-5 scripts/REFALT2haps.sh \
      --parfile ${PARFILE} --dir process/${PROJECT})
  echo "haps:    $jid_haps"

  # ── Step 5: smooth scan ───────────────────────────────────────────────────────
  jid_scan=$(sbatch --parsable --dependency=afterok:${jid_haps} \
      --array=1-5 scripts/haps2scan.freqsmooth.sh \
      --rfile ${DESIGN} --dir process/${PROJECT} --outdir ${SCAN} \
      --cov-smooth-kb ${COV_KB} --freq-smooth-kb ${FREQ_KB})
  echo "scan:    $jid_scan"

  # ── Step 6: concat + tarball ──────────────────────────────────────────────────
  jid_concat=$(sbatch --parsable --dependency=afterok:${jid_scan} \
      -A tdlong_lab -p standard --mem=10G \
      --wrap="bash scripts/concat_Chromosome_Scans.sh process/${PROJECT}/${SCAN}")
  echo "concat:  $jid_concat"

  # ── Step 7: publication figure ────────────────────────────────────────────────
  jid_fig=$(sbatch --parsable --dependency=afterok:${jid_concat} \
      -A tdlong_lab -p standard --mem=10G \
      --wrap="module load R/4.2.2 && Rscript ${FIGURE}")
  echo "figure:  $jid_fig"
