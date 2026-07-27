#!/bin/bash
#SBATCH --job-name=catalog_build
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00
# One CHROMOSOME per array task (--array=1-5). Submitted by build_catalog.sh.

# catalog_build.sh — annotate one chromosome's candidate founder SNPs.
#
# Part of the PROPOSED founder-catalog caller. A worker driven by build_catalog.sh.
# Calls the founders and emits EVERY candidate biallelic SNP with the annotations
# thresholds are decided on, so those thresholds become a fast downstream re-cut
# (catalog_filter.sh) instead of a full rebuild (XQTL2 #21):
#   CHROM POS REF ALT  dist_indel  AD_<founder1> AD_<founder2> ...
# where dist_indel = bp to the nearest founder indel, and AD_<f> = <ref>,<alt> for
# founder f. NO fixation / SnpGap filter here — catalog_filter.sh applies min-dp,
# maxaf, snpgap, and the segregation rule from these columns.
#
# Usage:
#   sbatch --array=1-5 catalog_build.sh --founders <founder_bams.txt> --catdir <catalog dir>

set -euo pipefail

REF=pipeline/ref/dm6.fa
THREADS=2

while [[ $# -gt 0 ]]; do
  case $1 in
    --founders) FOUNDER_BAMS="$2"; shift 2 ;;
    --catdir)   CATDIR="$2";       shift 2 ;;
    --ref)      REF="$2";          shift 2 ;;
    --threads)  THREADS="$2";      shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "${FOUNDER_BAMS:-}" ]] && { echo "Error: --founders required" >&2; exit 1; }
[[ -z "${CATDIR:-}" ]]      && { echo "Error: --catdir required" >&2; exit 1; }

chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$((SLURM_ARRAY_TASK_ID - 1))]:-}
[[ -z "$mychr" ]] && { echo "No chromosome for array index ${SLURM_ARRAY_TASK_ID:-<unset>}" >&2; exit 1; }

module load bcftools/1.21
mkdir -p "${CATDIR}/work"
raw="${CATDIR}/work/founders.calls.${mychr}.bcf"
annot="${CATDIR}/work/annot.${mychr}.txt"
indels="${CATDIR}/work/indels.${mychr}.txt"

# Call the founders on this chromosome, split multiallelics. This is DISCOVERY, so
# BAQ is ON (no -B): it suppresses false SNPs near indels, exactly what deciding
# whether a founder SNP is real should do (XQTL2 #23; matches the validated
# bam2bcf2REFALT.sh, which also leaves BAQ on). Counting samples is the opposite
# case and stays BAQ-off (catalog_count.sh, #22). No SnpGap filter here — near-indel
# SNPs are kept and annotated by distance instead.
bcftools mpileup -q 20 -Q 20 --max-depth 1000 -a FORMAT/AD,FORMAT/DP \
    -r "$mychr" -f "$REF" --threads "$THREADS" -b "$FOUNDER_BAMS" -Ou \
  | bcftools call -mv --threads "$THREADS" -Ou \
  | bcftools norm -f "$REF" -m - --threads "$THREADS" -Ob > "$raw"

# Founder indel positions on this chromosome (sorted), for the distance annotation.
bcftools view -v indels "$raw" | bcftools query -f '%POS\n' | sort -n > "$indels"

# Every candidate biallelic SNP + per-founder AD, annotated with distance to the
# nearest founder indel (two-pointer over the sorted indel list; input is
# coordinate-sorted so the pointer only advances).
bcftools view -m2 -M2 -v snps "$raw" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
  | awk -F'\t' -v OFS='\t' -v indelfile="$indels" '
    BEGIN { ni = 0; while ((getline p < indelfile) > 0) ind[++ni] = p + 0; j = 1 }
    {
        pos = $2 + 0
        while (j <= ni && ind[j] < pos) j++
        d = 1000000000
        if (j     <= ni) { dd = ind[j]   - pos; if (dd < d) d = dd }
        if (j - 1 >= 1)  { dd = pos - ind[j-1]; if (dd < d) d = dd }
        row = $1 OFS $2 OFS $3 OFS $4 OFS d
        for (i = 5; i <= NF; i++) row = row OFS $i
        print row
    }' > "$annot"

echo "${mychr}: $(wc -l < "$annot") candidate SNPs annotated -> $annot"
