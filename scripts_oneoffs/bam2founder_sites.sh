#!/bin/bash
#SBATCH --job-name=founderSites
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --array=1-5
#SBATCH --time=24:00:00

# Build a per-chromosome VCF site catalog from B-population founders only.
# Run once; the output is reused for all future B-pop projects.
#
# Usage: sbatch bam2founder_sites.sh <bam_list> <output_dir>
#   bam_list   : file listing founder BAM paths (e.g. helpfiles/B_founders.bams.txt)
#   output_dir : directory to write founder_sites.<chr>.vcf.gz  (e.g. process/B_founder_sites)
#
# Output per chromosome:
#   founder_sites.<chr>.vcf.gz      biallelic SNPs QUAL>59, founders only
#   founder_sites.<chr>.vcf.gz.tbi  tabix index

module load bcftools/1.21

ref="ref/dm6.fa"
bams=$1
output=$2

mkdir -p $output

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

bcftools mpileup -I -d 1000 -t $mychr -a "FORMAT/AD,FORMAT/DP" -f $ref -b $bams \
  | bcftools call -mv -Ob \
  | bcftools view -m2 -M2 -v snps -i 'QUAL>59' -Oz -o ${output}/founder_sites.${mychr}.vcf.gz

bcftools index -t ${output}/founder_sites.${mychr}.vcf.gz
