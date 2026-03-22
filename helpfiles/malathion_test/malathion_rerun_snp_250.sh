#!/bin/bash
# malathion_rerun_snp_250.sh — SNP scan only (250kb, smooth already exists)
#
# Run from XQTL2 project root:
#   bash helpfiles/malathion_test/malathion_rerun_snp_250.sh

bash scripts/run_snp_scan.sh \
    --design    helpfiles/malathion_test/design.txt \
    --dir       process/malathion_test \
    --scan      MALATHION_TEST_v2_smooth250 \
    --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founders  A1,A2,A3,A4,A5,A6,A7,AB8 \
    --figure    helpfiles/malathion_test/MALATHION_TEST_v2_smooth250_snp.R
