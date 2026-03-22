#!/bin/bash
# malathion_rerun_snp.sh
#
# Submit all 5 snp_scan chromosomes, then concat, then figures.
# Run from XQTL2 project root:
#   bash helpfiles/malathion_test/malathion_rerun_snp.sh

set -e

PROJECT=malathion_test
DESIGN=helpfiles/${PROJECT}/design.txt
SCAN=MALATHION_TEST_v2_smooth125
OUTDIR=process/${PROJECT}/${SCAN}
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8
SNP_TABLE=helpfiles/FREQ_SNPs_Apop.cM.txt.gz
FIGURE=helpfiles/${PROJECT}/MALATHION_TEST_v2_smooth125_snp.R

jid_snp=$(sbatch --parsable \
    --array=1-5 scripts/snp_scan.sh \
    --rfile     ${DESIGN} \
    --dir       process/${PROJECT} \
    --outdir    ${SCAN} \
    --snp-table ${SNP_TABLE} \
    --founders  ${FOUNDERS})
echo "snp_scan: $jid_snp"

jid_concat=$(sbatch --parsable \
    --dependency=afterok:${jid_snp} \
    -A tdlong_lab -p standard \
    --cpus-per-task=1 --mem-per-cpu=6G \
    --wrap="bash scripts/concat_snp_scans.sh ${OUTDIR}")
echo "concat:   $jid_concat"

jid_fig=$(sbatch --parsable \
    --dependency=afterok:${jid_concat} \
    -A tdlong_lab -p standard \
    --cpus-per-task=2 --mem-per-cpu=6G \
    --wrap="module load R/4.2.2 && Rscript ${FIGURE}")
echo "figure:   $jid_fig"

sbatch \
    --dependency=afterany:${jid_fig} \
    -A tdlong_lab -p standard \
    --cpus-per-task=1 --mem-per-cpu=6G \
    --wrap="$(cat <<SEFF
report=${OUTDIR}/snp_memory_report.txt
echo 'snp_scan memory report for ${SCAN}' > \$report
echo "Generated: \$(date)" >> \$report
for t in 1 2 3 4 5; do
    echo "" >> \$report
    echo "=== snp_scan ${jid_snp}_\$t ===" >> \$report
    seff ${jid_snp}_\$t >> \$report
done
echo "" >> \$report
echo "=== concat ${jid_concat} ===" >> \$report
seff ${jid_concat} >> \$report
echo "" >> \$report
echo "=== figure ${jid_fig} ===" >> \$report
seff ${jid_fig} >> \$report
SEFF
)"
echo "seff report -> ${OUTDIR}/snp_memory_report.txt"
