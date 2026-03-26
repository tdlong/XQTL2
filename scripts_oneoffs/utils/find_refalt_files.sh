#!/bin/bash
# find_refalt_files.sh
# Run from anywhere on the cluster
# bash scripts/find_refalt_files.sh 2>&1 | tee refalt_report.txt

SEP="════════════════════════════════════════════════════════════════"

echo "$SEP"
echo "MALATHION — RefAlt files in newpipeline_Nov23"
echo "$SEP"
find /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23 -name "RefAlt*" 2>/dev/null | sort
echo ""

echo "$SEP"
echo "MALATHION — process subdirs in newpipeline_Nov23"
echo "$SEP"
ls /dfs7/adl/tdlong/fly_pool/newpipeline_Nov23/ 2>/dev/null
echo ""

echo "$SEP"
echo "PUPATION — RefAlt files in Sarah's pupal_dynamic"
echo "$SEP"
find /dfs7/adl/sruckman/XQTL/XQTL2/process/pupal_dynamic -name "RefAlt*" 2>/dev/null | sort
echo ""

echo "$SEP"
echo "PUPATION — hap files (R.haps.*) in Sarah's pupal_dynamic"
echo "$SEP"
find /dfs7/adl/sruckman/XQTL/XQTL2/process/pupal_dynamic -name "R.haps.*" 2>/dev/null | sort
echo ""

echo "$SEP"
echo "AGING — hap files in process/AGE_SY"
echo "$SEP"
ls /dfs7/adl/tdlong/fly_pool/XQTL2/process/AGE_SY/*.rds 2>/dev/null | sort
echo ""

echo "$SEP"
echo "AGING — scan job status (AGE_SY20_M_freqs250)"
echo "$SEP"
ls /dfs7/adl/tdlong/fly_pool/XQTL2/process/AGE_SY/AGE_SY20_M_freqs250/ 2>/dev/null | head -20
echo ""
