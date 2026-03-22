#!/bin/bash
# malathion_test_v2_smooth250.sh
#
# Run hap scan pipeline with 250 kb smoothing.
# Starts from smooth_haps (raw haps RDS already exist).
# Run from XQTL2 project root:
#   bash helpfiles/malathion_test/malathion_test_v2_smooth250.sh

set -e

PROJECT=malathion_test
DESIGN=helpfiles/${PROJECT}/design.txt
SCAN=MALATHION_TEST_v2_smooth250
OUTDIR=process/${PROJECT}/${SCAN}
FIGURE=helpfiles/${PROJECT}/MALATHION_TEST_v2_smooth250.R

mkdir -p ${OUTDIR}

jid_smooth=$(sbatch --parsable \
    --array=1-5 scripts_freqsmooth/smooth_haps.sh \
    --rfile     ${DESIGN} \
    --dir       process/${PROJECT} \
    --outdir    ${SCAN} \
    --smooth-kb 250)
echo "smooth: $jid_smooth"

jid_scan=$(sbatch --parsable \
    --dependency=afterok:${jid_smooth} \
    --array=1-5 scripts_freqsmooth/freqsmooth_scan.sh \
    --rfile  ${DESIGN} \
    --dir    process/${PROJECT} \
    --outdir ${SCAN})
echo "scan:   $jid_scan"

jid_concat=$(sbatch --parsable \
    --dependency=afterok:${jid_scan} \
    -A tdlong_lab -p standard \
    --cpus-per-task=2 --mem-per-cpu=6G \
    --wrap="bash scripts/concat_Chromosome_Scans.sh ${OUTDIR}")
echo "concat: $jid_concat"

jid_fig=$(sbatch --parsable \
    --dependency=afterok:${jid_concat} \
    -A tdlong_lab -p standard \
    --cpus-per-task=2 --mem-per-cpu=6G \
    --wrap="module load R/4.2.2 && Rscript ${FIGURE}")
echo "figure: $jid_fig"

sbatch \
    --dependency=afterany:${jid_fig} \
    -A tdlong_lab -p standard \
    --cpus-per-task=1 --mem-per-cpu=6G \
    --wrap="$(cat <<SEFF
report=${OUTDIR}/memory_report.txt
echo 'Memory report for ${SCAN}' > \$report
echo "Generated: \$(date)" >> \$report
for t in 1 2 3 4 5; do
    echo "" >> \$report
    echo "=== smooth ${jid_smooth}_\$t ===" >> \$report
    seff ${jid_smooth}_\$t >> \$report
    echo "=== scan ${jid_scan}_\$t ===" >> \$report
    seff ${jid_scan}_\$t >> \$report
done
echo "" >> \$report
echo "=== concat ${jid_concat} ===" >> \$report
seff ${jid_concat} >> \$report
echo "" >> \$report
echo "=== figure ${jid_fig} ===" >> \$report
seff ${jid_fig} >> \$report
SEFF
)"
echo "seff report -> ${OUTDIR}/memory_report.txt"
