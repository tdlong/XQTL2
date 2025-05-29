#!/bin/bash
#SBATCH --job-name=callSNPs
#SBATCH -A tdlong_lab        ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --cpus-per-task=2 
#SBATCH --array=1-5
#SBATCH --time=5-00:00:00

module load bcftools/1.21

ref="ref/dm6.fa"
# passed from command line
bams=$1
output=$2

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}
# assume I can write to process
bcftools mpileup -I -d 1000 -t $mychr -a "FORMAT/AD,FORMAT/DP" -f $ref -b $bams | bcftools call -mv -Ob > ${output}/calls.$mychr.bcf  
echo -ne "CHROM\tPOS" > ${output}/RefAlt.$mychr.txt
bcftools query -l ${output}/calls.$mychr.bcf | awk '{printf("\tREF_%s\tALT_%s",$1,$1)}' >> ${output}/RefAlt.$mychr.txt
echo -ne "\n" >> ${output}/RefAlt.$mychr.txt
bcftools view -m2 -M2 -v snps -i 'QUAL>59' ${output}/calls.$mychr.bcf | bcftools query -f'%CHROM %POS [ %AD{0} %AD{1}] [%GT]\n' | grep -v '\.' | awk 'NF-=1' >>${output}/RefAlt.$mychr.txt


