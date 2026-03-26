#!/bin/bash
#SBATCH --job-name=merge_b3852
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1

module load samtools/1.15.1

samtools merge data/founders/b3852.merged.bam \
  /dfs7/adl/tdlong/check_fly_strains/data/bam/June_2025/b3852.bam \
  /dfs7/adl/tdlong/check_fly_strains/data/bam/June_2025/b3852_r.bam

samtools index data/founders/b3852.merged.bam
