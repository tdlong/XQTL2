#!/bin/bash
# malathion_test_v2.sh — full scan pipeline for the malathion test dataset
#
# Run from XQTL2 project root:
#   bash helpfiles/malathion_test/malathion_test_v2.sh

bash scripts/run_scan.sh \
    --design    helpfiles/malathion_test/design.txt \
    --dir       process/malathion_test \
    --scan      MALATHION_TEST_v2_smooth250 \
    --smooth    250 \
    --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founders  A1,A2,A3,A4,A5,A6,A7,AB8 \
    --figure    helpfiles/malathion_test/MALATHION_TEST_v2_smooth250.R
