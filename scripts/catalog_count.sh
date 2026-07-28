#!/bin/bash
#SBATCH --job-name=catalog_count
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-00:00:00
# One SAMPLE per array task (--array=1-<#BAMs>). Submitted by call_samples.sh.
# 1 core: mpileup is single-threaded and uses ~1 core / <300 MB.

# catalog_count.sh — count REF/ALT reads for ONE BAM against a catalog.
#
# Part of the PROPOSED founder-catalog caller; a worker driven by call_samples.sh.
# Counts the BAM at line $SLURM_ARRAY_TASK_ID of <bam_list> against the catalog in
# <catalog dir>, writing <calls dir>/counts/<SM>.tsv.gz. No skip guard: it counts
# exactly the BAM it is given and overwrites. (To ADD samples, call_samples.sh is
# run with a bam list of just the new BAMs — their counts land alongside the rest.)
#
# Counting is deterministic (BAQ off -B; positions fixed by -T catalog), so counts
# do not depend on interval spec or on which other samples are in the run. These
# are POOLED samples: we read REF/ALT depths straight from mpileup and do NOT run a
# genotype caller (see below / XQTL2 #22). Column name = the BAM's SM read-group tag.
#
# Usage:
#   sbatch --array=1-<#BAMs> catalog_count.sh <bam_list> <catalog dir> <calls dir>

set -euo pipefail

ref="pipeline/ref/dm6.fa"
bams=$1
catdir=$2
callsdir=$3
cat="${catdir}/catalog.tsv.gz"

[[ -f "$cat" ]] || { echo "Error: catalog not found: $cat (build it first)" >&2; exit 1; }

# NB: bcftools/1.21 and samtools/1.10 link incompatible htslib versions and must
# never be loaded in the same shell. Load bcftools here; read the SM tag with
# samtools in an isolated subshell below.
module load bcftools/1.21
mkdir -p "${callsdir}/counts"

bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$bams")
[[ -z "$bam" ]] && { echo "No BAM at line $SLURM_ARRAY_TASK_ID of $bams" >&2; exit 1; }

name=$( module load samtools/1.10; samtools view -H "$bam" | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print substr($i,4)}' | sort -u | head -1 )
[[ -z "$name" ]] && { echo "No SM tag in $bam" >&2; exit 1; }
out="${callsdir}/counts/${name}.tsv.gz"

# These are POOLED samples, not diploid individuals — the task is to COUNT REF vs
# ALT reads at a known catalog SNP, not to genotype it. Take AD straight from
# `bcftools mpileup`. Do NOT route it through `bcftools call`: the diploid genotype
# model calls a per-sample-rare minor allele "homozygous REF" and erases its real
# ALT reads (XQTL2 #22). mpileup's own -q/-Q/BAQ already exclude bad *reads*, which
# is the only legitimate exclusion here.
#
# AD from mpileup is [REF, <first ALT>, <*>]; for a clean biallelic catalog SNP the
# first ALT is the catalog ALT (or 0 when the sample has no ALT reads). Merge keys
# on (CHROM,POS) only, so no genotype step is needed to keep alleles consistent.
#
# -I (--skip-indels): without it, at a catalog SNP where the sample also carries an
# indel, mpileup emits BOTH a SNP record and an indel record at the same POS, so the
# count file has duplicate positions and the merge cartesian-explodes (XQTL2 #27).
# -I suppresses indel *calling* only; the reads piled at the SNP (its REF/ALT depths)
# are unchanged. The validated bam2bcf2REFALT.sh does the same.
{
  printf 'CHROM\tPOS\tREF_%s\tALT_%s\n' "$name" "$name"
  bcftools mpileup -I -B -q 20 -Q 20 --max-depth 2000 -T "$cat" -a FORMAT/AD \
      -f "$ref" "$bam" -Ou \
    | bcftools query -f '%CHROM\t%POS[\t%AD]\n' \
    | awk -F'\t' -v OFS='\t' '{ split($3,a,","); print $1,$2,(a[1]==""||a[1]=="."?0:a[1]),(a[2]==""||a[2]=="."?0:a[2]) }'
} | bgzip > "$out"

echo "counted ${name} -> ${out}"
