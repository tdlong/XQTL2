#!/bin/bash
#SBATCH --job-name=bwa_B6885
#SBATCH -A tdlong_lab        ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --cpus-per-task=4

module load bwa/0.7.17
module load samtools/1.10
module load bcftools/1.21
module load java/17
module load picard-tools/1.87

# Reference genome
ref="ref/dm6.fa"

# Arguments:
#   $1 = mapping file (tab-separated: fastq_prefix  shortname)
#   $2 = path to raw FASTQ directory
#   $3 = output BAM directory
files=$1
INDIR=$2
OUTDIR=$3

# Read the mapping file for this array task
prefix=$(head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f1)
shortname=$(head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f2)

# Illumina-style FASTQ naming: {prefix}_R1_001.fastq.gz
R1="${INDIR}/${prefix}_R1_001.fastq.gz"
R2="${INDIR}/${prefix}_R2_001.fastq.gz"

echo "Aligning $shortname: $R1 + $R2"

# Align, sort, add read groups, index (identical to fq2bam.sh)
bwa mem -t 4 -M $ref ${R1} ${R2} | samtools view -bS - > ${OUTDIR}/$shortname.temp1.bam
samtools sort ${OUTDIR}/$shortname.temp1.bam -o ${OUTDIR}/$shortname.temp2.bam
rm ${OUTDIR}/$shortname.temp1.bam

prog="/opt/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar"
java -Xmx20g -jar $prog \
    I=${OUTDIR}/$shortname.temp2.bam \
    O=${OUTDIR}/$shortname.bam \
    SORT_ORDER=coordinate \
    RGPL=illumina \
    RGPU=D109LACXX \
    RGLB=Lib1 \
    RGID=$shortname \
    RGSM=$shortname \
    VALIDATION_STRINGENCY=LENIENT

samtools index ${OUTDIR}/$shortname.bam
rm ${OUTDIR}/$shortname.temp2.bam

echo "Done: ${OUTDIR}/$shortname.bam"
