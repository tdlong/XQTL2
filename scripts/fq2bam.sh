#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -A tdlong_lab        ## account to charge 
#SBATCH -p standard          ## partition/queue name
#SBATCH --cpus-per-task=4 

module load bwa/0.7.17
module load samtools/1.10
module load bcftools/1.21
module load java/17 
module load picard-tools/1.87

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
R1=`ls $INDIR/*.txt.gz | grep READ1 | grep $BC`
R2=`echo $R1 | sed 's/READ1/READ2/'`

bwa mem -t 4 -M $ref ${R1} ${R2} | samtools view -bS - > ${OUTDIR}/$shortname.temp1.bam
samtools sort ${OUTDIR}/$shortname.temp1.bam -o ${OUTDIR}/$shortname.temp2.bam
rm ${OUTDIR}/$shortname.temp1.bam
# Locate AddOrReplaceReadGroups.jar via the module environment.
# After "module load picard-tools", common variable names are $PICARD (full jar path)
# or $PICARD_HOME / $PICARDPATH (installation directory).  Try each in order.
if   [[ -f "$PICARD" ]];                                            then prog="$PICARD"
elif [[ -n "$PICARD_HOME"  ]] && find "$PICARD_HOME"  -name "AddOrReplaceReadGroups.jar" -quit 2>/dev/null; then
     prog=$(find "$PICARD_HOME"  -name "AddOrReplaceReadGroups.jar" | head -1)
elif [[ -n "$PICARDPATH"   ]] && find "$PICARDPATH"   -name "AddOrReplaceReadGroups.jar" -quit 2>/dev/null; then
     prog=$(find "$PICARDPATH"   -name "AddOrReplaceReadGroups.jar" | head -1)
else
     echo "Error: cannot locate AddOrReplaceReadGroups.jar."
     echo "       Load the picard-tools module and verify one of these is set:"
     echo "         \$PICARD, \$PICARD_HOME, \$PICARDPATH"
     exit 1
fi
java -Xmx20g -jar $prog  I=${OUTDIR}/$shortname.temp2.bam O=${OUTDIR}/$shortname.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=$shortname RGSM=$shortname VALIDATION_STRINGENCY=LENIENT
samtools index ${OUTDIR}/$shortname.bam
rm ${OUTDIR}/$shortname.temp2.bam

