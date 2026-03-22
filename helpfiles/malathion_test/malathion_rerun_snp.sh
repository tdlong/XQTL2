#!/bin/bash
# malathion_rerun_snp.sh — SNP scan only (125kb, smooth already exists)
#
# Run from XQTL2 project root:
#   bash helpfiles/malathion_test/malathion_rerun_snp.sh

bash scripts/run_scan.sh \
    --design    helpfiles/malathion_test/design.txt \
    --dir       process/malathion_test \
    --scan      MALATHION_TEST_v2_smooth125 \
    --smooth    125 \
    --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founders  A1,A2,A3,A4,A5,A6,A7,AB8 \
    --figure    helpfiles/malathion_test/MALATHION_TEST_v2_smooth125_snp.R
