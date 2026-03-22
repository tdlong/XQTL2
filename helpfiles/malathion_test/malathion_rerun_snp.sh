#!/bin/bash
# malathion_rerun_snp.sh
#
# Cancel any running snp_scan jobs and resubmit all 5 chromosomes.
# Run from XQTL2 project root:
#   bash helpfiles/malathion_test/malathion_rerun_snp.sh

set -e

PROJECT=malathion_test
DESIGN=helpfiles/${PROJECT}/design.txt
SCAN=MALATHION_TEST_v2_smooth125
OUTDIR=process/${PROJECT}/${SCAN}
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8
SNP_TABLE=helpfiles/FREQ_SNPs_Apop.cM.txt.gz

jid=$(sbatch --parsable \
    --array=1-5 scripts_freqsmooth/snp_scan.sh \
    --rfile     ${DESIGN} \
    --dir       process/${PROJECT} \
    --outdir    ${SCAN} \
    --snp-table ${SNP_TABLE} \
    --founders  ${FOUNDERS})
echo "snp_scan: $jid"

# seff report after all tasks finish
sbatch \
    --dependency=afterany:${jid} \
    -A tdlong_lab -p standard \
    --cpus-per-task=1 --mem-per-cpu=6G \
    --wrap="$(cat <<SEFF
report=${OUTDIR}/snp_memory_report.txt
echo 'snp_scan memory report' > \$report
echo "Generated: \$(date)" >> \$report
for t in 1 2 3 4 5; do
    echo "" >> \$report
    echo "=== snp_scan ${jid}_\$t ===" >> \$report
    seff ${jid}_\$t >> \$report
done
SEFF
)"
echo "seff report -> ${OUTDIR}/snp_memory_report.txt"
