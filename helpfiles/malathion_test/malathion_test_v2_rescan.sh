#!/bin/bash
# malathion_test_v2_rescan.sh — rerun hap scan from smooth onward (125kb)
#
# Run from XQTL2 project root:
#   bash helpfiles/malathion_test/malathion_test_v2_rescan.sh

bash scripts/run_scan.sh \
    --design helpfiles/malathion_test/design.txt \
    --dir    process/malathion_test \
    --scan   MALATHION_TEST_v2_smooth125 \
    --smooth 125

# Figures (run separately after scan completes):
#   Rscript helpfiles/malathion_test/MALATHION_TEST_v2_smooth125.R
