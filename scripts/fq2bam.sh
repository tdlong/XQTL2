#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -A tdlong_lab        ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-00:00:00

module load bwa/0.7.17
module load samtools/1.10
module load bcftools/1.21
module load picard-tools/3.3.0

# assume this exists
ref="ref/dm6.fa"
# from command line = Barcode A, Barcode B, sample name
files=$1
# this is the path to the raw data
INDIR=$2
# this is where we will write the bams to
OUTDIR=$3


BCA=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f1`
BCB=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f2`
shortname=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f3`
BC="${BCA}-${BCB}"
# auto-detect naming convention: old (*.txt.gz, READ1/READ2) or new (*.fastq.gz, -R1/-R2)
R1=`ls $INDIR/*.txt.gz 2>/dev/null | grep READ1 | grep $BC` || true
if [[ -n "$R1" ]]; then
    R2=`echo $R1 | sed 's/READ1/READ2/'`
else
    R1=`ls $INDIR/*.fastq.gz 2>/dev/null | grep "\-R1\." | grep $BC` || true
    R2=`echo $R1 | sed 's/-R1\./-R2./'`
fi
if [[ -z "$R1" ]]; then
    echo "Error: no R1 fastq found for barcode $BC in $INDIR" >&2
    exit 1
fi

bwa mem -t 4 -M $ref ${R1} ${R2} | samtools view -bS - > ${OUTDIR}/$shortname.temp1.bam
samtools sort ${OUTDIR}/$shortname.temp1.bam -o ${OUTDIR}/$shortname.temp2.bam
rm ${OUTDIR}/$shortname.temp1.bam
java -Xmx20g -jar $PICARDTOOLS/picard.jar AddOrReplaceReadGroups I=${OUTDIR}/$shortname.temp2.bam O=${OUTDIR}/$shortname.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=$shortname RGSM=$shortname VALIDATION_STRINGENCY=LENIENT
samtools index ${OUTDIR}/$shortname.bam
rm ${OUTDIR}/$shortname.temp2.bam

