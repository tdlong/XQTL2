#!/bin/bash
###############################################################################
# chrX_check.sh — Rerun chrX scan steps for ZINC_Hanson_v3 and ZINC2_F_v3,
#                 then generate a chrX-only SNP Manhattan figure.
#
# Chain: smooth (2x6G) → hap_scan → snp_scan → figure
# Output: output/chrX_snp_manhattan.png
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

set -e
WD=/dfs7/adl/tdlong/fly_pool/XQTL2
mkdir -p output logs

# ── smooth chrX (2x6G — chrX needs 12G) ─────────────────────────────────────
jid_smooth_pa=$(sbatch --parsable \
    -A tdlong_lab -p standard --cpus-per-task=2 --mem-per-cpu=6G --time=1:00:00 \
    --job-name=smooth_chrX_pa \
    --output=logs/smooth_chrX_pa.out \
    --wrap="module load R/4.2.2 && Rscript scripts/smooth_haps.R \
        --chr chrX --dir process/ZINC_Hanson --outdir ZINC_Hanson_v3 \
        --rfile helpfiles/ZINC_Hanson/ZINC_Hanson.test.txt --smooth-kb 250")
echo "smooth_chrX_pa:   ${jid_smooth_pa}"

jid_smooth_pb=$(sbatch --parsable \
    -A tdlong_lab -p standard --cpus-per-task=2 --mem-per-cpu=6G --time=1:00:00 \
    --job-name=smooth_chrX_pb \
    --output=logs/smooth_chrX_pb.out \
    --wrap="module load R/4.2.2 && Rscript scripts/smooth_haps.R \
        --chr chrX --dir process/ZINC2 --outdir ZINC2_F_v3 \
        --rfile helpfiles/ZINC2/Zinc2.test.F.txt --smooth-kb 250")
echo "smooth_chrX_pb:   ${jid_smooth_pb}"

# ── hap_scan chrX ─────────────────────────────────────────────────────────────
jid_hap_pa=$(sbatch --parsable \
    --dependency=afterok:${jid_smooth_pa} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
    --job-name=hap_chrX_pa \
    --output=logs/hap_chrX_pa.out \
    --wrap="module load R/4.2.2 && Rscript scripts/hap_scan.R \
        --chr chrX --dir process/ZINC_Hanson/ZINC_Hanson_v3 \
        --outdir ZINC_Hanson_v3 --rfile helpfiles/ZINC_Hanson/ZINC_Hanson.test.txt")
echo "hap_chrX_pa:      ${jid_hap_pa}"

jid_hap_pb=$(sbatch --parsable \
    --dependency=afterok:${jid_smooth_pb} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
    --job-name=hap_chrX_pb \
    --output=logs/hap_chrX_pb.out \
    --wrap="module load R/4.2.2 && Rscript scripts/hap_scan.R \
        --chr chrX --dir process/ZINC2/ZINC2_F_v3 \
        --outdir ZINC2_F_v3 --rfile helpfiles/ZINC2/Zinc2.test.F.txt")
echo "hap_chrX_pb:      ${jid_hap_pb}"

# ── snp_scan chrX ─────────────────────────────────────────────────────────────
jid_snp_pa=$(sbatch --parsable \
    --dependency=afterok:${jid_hap_pa} \
    -A tdlong_lab -p standard --cpus-per-task=2 --mem-per-cpu=6G --time=1:00:00 \
    --job-name=snp_chrX_pa \
    --output=logs/snp_chrX_pa.out \
    --wrap="module load R/4.2.2 && Rscript scripts/snp_scan.R \
        --chr chrX --dir process/ZINC_Hanson/ZINC_Hanson_v3 \
        --outdir ZINC_Hanson_v3 --rfile helpfiles/ZINC_Hanson/ZINC_Hanson.test.txt \
        --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
        --founders A1,A2,A3,A4,A5,A6,A7,AB8")
echo "snp_chrX_pa:      ${jid_snp_pa}"

jid_snp_pb=$(sbatch --parsable \
    --dependency=afterok:${jid_hap_pb} \
    -A tdlong_lab -p standard --cpus-per-task=2 --mem-per-cpu=6G --time=1:00:00 \
    --job-name=snp_chrX_pb \
    --output=logs/snp_chrX_pb.out \
    --wrap="module load R/4.2.2 && Rscript scripts/snp_scan.R \
        --chr chrX --dir process/ZINC2/ZINC2_F_v3 \
        --outdir ZINC2_F_v3 --rfile helpfiles/ZINC2/Zinc2.test.F.txt \
        --snp-table helpfiles/FREQ_SNPs_Bpop.cM.txt.gz \
        --founders B1,B2,B3,B4,B5,B6,B7,AB8")
echo "snp_chrX_pb:      ${jid_snp_pb}"

# ── chrX SNP Manhattan figure ─────────────────────────────────────────────────
jid_fig=$(sbatch --parsable \
    --dependency=afterok:${jid_snp_pa},afterok:${jid_snp_pb} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=0:30:00 \
    --job-name=chrX_fig \
    --output=logs/chrX_fig.out \
    --wrap="module load R/4.2.2 && cd ${WD} && Rscript scripts_oneoffs/ZINC2/chrX_snp_manhattan.R")
echo "chrX_fig:         ${jid_fig}"

echo ""
echo "When complete, download:"
echo "  scp tdlong@hpc3.rcic.uci.edu:${WD}/output/chrX_snp_manhattan.png ."
