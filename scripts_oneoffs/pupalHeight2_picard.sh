#!/bin/bash
#SBATCH --job-name=picard
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2:00:00

module load picard-tools/3.3.0
module load samtools/1.10

files=helpfiles/pupalHeight2/pupalHeight2.barcodes.txt
OUTDIR=data/bam/pupalHeight2

shortname=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f3`

java -Xmx20g -jar $PICARDTOOLS/picard.jar AddOrReplaceReadGroups I=${OUTDIR}/$shortname.temp2.bam O=${OUTDIR}/$shortname.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=$shortname RGSM=$shortname VALIDATION_STRINGENCY=LENIENT
samtools index ${OUTDIR}/$shortname.bam
rm ${OUTDIR}/$shortname.temp2.bam
