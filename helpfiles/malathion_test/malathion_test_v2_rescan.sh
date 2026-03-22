#!/bin/bash
# malathion_test_v2_rescan.sh
#
# Rerun steps 5b onward (smooth output already exists).
# Run from XQTL2 project root:
#   bash helpfiles/malathion_test/malathion_test_v2_rescan.sh

set -e

PROJECT=malathion_test
DESIGN=helpfiles/${PROJECT}/design.txt
SCAN=MALATHION_TEST_v2_smooth125
FIGURE=helpfiles/${PROJECT}/MALATHION_TEST_v2_smooth125.R

jid_scan=$(sbatch --parsable \
    --array=1-5 scripts_freqsmooth/freqsmooth_scan.sh \
    --rfile  ${DESIGN} \
    --dir    process/${PROJECT} \
    --outdir ${SCAN})
echo "scan:   $jid_scan"

jid_concat=$(sbatch --parsable \
    --dependency=afterok:${jid_scan} \
    -A tdlong_lab -p standard \
    --cpus-per-task=2 --mem-per-cpu=6G \
    --wrap="bash scripts/concat_Chromosome_Scans.sh process/${PROJECT}/${SCAN}")
echo "concat: $jid_concat"

jid_fig=$(sbatch --parsable \
    --dependency=afterok:${jid_concat} \
    -A tdlong_lab -p standard \
    --cpus-per-task=2 --mem-per-cpu=6G \
    --wrap="module load R/4.2.2 && Rscript ${FIGURE}")
echo "figure: $jid_fig"
