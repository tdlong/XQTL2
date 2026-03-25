#!/bin/bash
#SBATCH --job-name=merge_bams
#SBATCH -A tdlong_lab        ## account to charge 
#SBATCH -p standard          ## partition/queue name

module load samtools/1.10
module load java/17 
module load picard-tools/1.87

# Define directories
TEMP_A="data/bam/AGE_SY_tempA"
TEMP_B="data/bam/AGE_SY_tempB"
OUTPUT_DIR="data/bam/AGE_SY"

echo "Starting BAM merge process at $(date)"

# Loop through all BAM files in tempA
for BAM_A in ${TEMP_A}/*.bam; do
    # Get sample name
    SAMPLE=$(basename ${BAM_A} .bam)
    
    BAM_B="${TEMP_B}/${SAMPLE}.bam"
    OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE}.bam"
    
    echo "----------------------------------------"
    echo "Processing: ${SAMPLE}"
    
    # Check if corresponding file exists in tempB
    if [[ ! -f ${BAM_B} ]]; then
        echo "WARNING: ${BAM_B} not found, skipping ${SAMPLE}"
        continue
    fi
    
    # Merge using samtools (faster for this loop approach)
    echo "Merging BAM files..."
    samtools merge -f ${OUTPUT_BAM} ${BAM_A} ${BAM_B}
    
    # Index the merged BAM
    echo "Indexing merged BAM..."
    samtools index ${OUTPUT_BAM}
    
    # Quick validation
    if [[ -f ${OUTPUT_BAM} && -f ${OUTPUT_BAM}.bai ]]; then
        echo "✓ Successfully merged and indexed ${SAMPLE}"
        samtools flagstat ${OUTPUT_BAM} | head -n 1
    else
        echo "✗ ERROR: Failed to create output for ${SAMPLE}"
    fi
done

