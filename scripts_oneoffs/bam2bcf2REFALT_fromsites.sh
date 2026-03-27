#!/bin/bash
#SBATCH --job-name=refaltFromSites
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-5
#SBATCH --time=12:00:00

# Generate a per-sample RefAlt file using a pre-built founder site catalog.
# Run once per sample (not once per project). To add new samples to an existing
# project, submit this script for the new samples only, then re-run consolidation.
#
# Usage: sbatch bam2bcf2REFALT_fromsites.sh <bam> <sites_dir> <output_dir>
#   bam        : single sample BAM file
#   sites_dir  : directory containing founder_sites.<chr>.vcf.gz (from bam2founder_sites.sh)
#   output_dir : directory to write RefAlt.<samplename>.<chr>.txt
#
# Output per chromosome:
#   RefAlt.<samplename>.<chr>.txt   tab-delimited: CHROM POS REF_<name> ALT_<name>
#
# NOTE: No QUAL filter applied here — sites are pre-validated from founders.
#       Sites with zero alt reads in this sample will be absent (handled by
#       consolidate_refalt_fromsites.R which fills missing sites with 0).

module load bcftools/1.21

ref="ref/dm6.fa"
bam=$1
sites_dir=$2
output=$3

mkdir -p $output

# Sample name from BAM header SM tag via a minimal mpileup — avoids loading samtools
# alongside bcftools (htslib version conflict). bcftools always emits the VCF header
# with the SM tag in the #CHROM line even when the region has no reads.
samplename=$(bcftools mpileup -f $ref -r chrX:1-100 $bam 2>/dev/null | awk '/^#CHROM/{print $NF; exit}')

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

sites=${sites_dir}/founder_sites.${mychr}.vcf.gz
sites_tab=${sites_dir}/founder_sites.${mychr}.tab

outfile=${output}/RefAlt.${samplename}.${mychr}.txt

# Write header
echo -e "CHROM\tPOS\tREF_${samplename}\tALT_${samplename}" > $outfile

# Pileup restricted to founder sites only, then call variants in this sample
bcftools mpileup -I -d 1000 -t $mychr -T $sites -a "FORMAT/AD,FORMAT/DP" -f $ref $bam \
  | bcftools call -C alleles -T $sites_tab -m -Ob \
  | bcftools query -f'%CHROM\t%POS\t[%AD{0}\t%AD{1}]\n' \
  | grep -v '\.' \
  >> $outfile
