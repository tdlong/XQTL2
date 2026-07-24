#!/bin/bash
#SBATCH --job-name=catalog_build
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00
# One CHROMOSOME per array task (--array=1-5). Submitted by build_catalog.sh.

# catalog_build.sh — build one chromosome's catalog piece from the founder BAMs.
#
# Part of the PROPOSED founder-catalog caller (separate from the validated
# bam2bcf2REFALT.sh). A worker, normally driven by build_catalog.sh. Writes into
# <catdir>/work/. No reuse guard: running it (re)builds the piece.
#
# A biallelic SNP is KEPT if, across the founders (an --exempt founder on this
# chromosome is skipped entirely — as if not a founder):
#   1. every founder has depth >= --min-dp,
#   2. every founder is near-fixed (ALT freq <= --maxaf or >= 1-maxaf),
#   3. it segregates (at least one founder REF-fixed and one ALT-fixed),
#   4. it is not within --snpgap bp of a founder indel (bcftools --SnpGap).
# A per-rule tally of how many candidates fail each rule is written to
# <catdir>/work/stats.<chr>.txt (build_catalog.sh sums it).
#
# Usage:
#   sbatch --array=1-5 catalog_build.sh --founders <founder_bams.txt> --catdir <catalog dir> \
#          [--min-dp 10 --maxaf 0.03 --snpgap 5 --exempt-founders B5:chr2L]

set -euo pipefail

REF=pipeline/ref/dm6.fa
MIN_DP=10
MAXAF=0.03
SNPGAP=5
THREADS=2
EXEMPT="B5:chr2L"

while [[ $# -gt 0 ]]; do
  case $1 in
    --founders)        FOUNDER_BAMS="$2"; shift 2 ;;
    --catdir)          CATDIR="$2";       shift 2 ;;
    --ref)             REF="$2";          shift 2 ;;
    --min-dp)          MIN_DP="$2";       shift 2 ;;
    --maxaf)           MAXAF="$2";        shift 2 ;;
    --snpgap)          SNPGAP="$2";       shift 2 ;;
    --exempt-founders) EXEMPT="$2";       shift 2 ;;
    --threads)         THREADS="$2";      shift 2 ;;
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
bed="${CATDIR}/work/catalog.${mychr}.bed"
stats="${CATDIR}/work/stats.${mychr}.txt"

# Call the founders on this chromosome (BAQ off), split multiallelics, flag SNPs near indels.
bcftools mpileup -B -q 20 -Q 20 --max-depth 1000 -a FORMAT/AD,FORMAT/DP \
    -r "$mychr" -f "$REF" --threads "$THREADS" -b "$FOUNDER_BAMS" -Ou \
  | bcftools call -mv --threads "$THREADS" -Ou \
  | bcftools norm -f "$REF" -m - --threads "$THREADS" -Ou \
  | bcftools filter -g "$SNPGAP" --threads "$THREADS" -Ob > "$raw"

# Founder order == the AD column order below. Turn --exempt-founders (NAME or
# NAME:CHR) into the AD fields to skip on THIS chromosome.
fnames=( $(bcftools query -l "$raw") )
skipcols=""
for i in "${!fnames[@]}"; do
  for e in ${EXEMPT//,/ }; do
    ename="${e%%:*}"; echr=""; [[ "$e" == *:* ]] && echr="${e#*:}"
    if [[ "${fnames[$i]}" == "$ename" ]] && { [[ -z "$echr" ]] || [[ "$echr" == "$mychr" ]]; }; then
      skipcols="${skipcols} $((5 + i))"
    fi
  done
done
[[ -n "${skipcols// }" ]] && echo "${mychr}: exempting founder field(s)${skipcols} (${EXEMPT})"

# Fixation filter -> bed (CHROM POS0 POS REF ALT); per-rule tally -> stats file.
bcftools view -m2 -M2 -v snps -e 'FILTER~"SnpGap"' "$raw" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
  | awk -F'\t' -v OFS='\t' -v mindp="$MIN_DP" -v maxaf="$MAXAF" \
        -v skipcols="$skipcols" -v statsfile="$stats" '
    BEGIN { n = split(skipcols, s, " "); for (k = 1; k <= n; k++) if (s[k] != "") skip[s[k]] = 1 }
    {
        total++
        segR = 0; segA = 0; reason = ""
        for (i = 5; i <= NF; i++) {
            if (i in skip) continue
            split($i, a, ",")
            r = a[1] + 0; v = a[2] + 0; dp = r + v
            if (dp < mindp) { reason = "depth"; break }
            af = v / dp
            if (af > maxaf && af < 1 - maxaf) { reason = "notfixed"; break }
            if (af <= maxaf)     segR = 1
            if (af >= 1 - maxaf) segA = 1
        }
        if (reason == "depth")      { c_depth++;    next }
        if (reason == "notfixed")   { c_notfixed++; next }
        if (!(segR && segA))        { c_notseg++;   next }
        c_kept++
        print $1, $2 - 1, $2, $3, $4
    }
    END {
        printf "candidates=%d depth_fail=%d notfixed_fail=%d notseg_fail=%d kept=%d\n",
               total, c_depth, c_notfixed, c_notseg, c_kept > statsfile
    }' > "$bed"

echo "${mychr}: $(cat "$stats")  -> $bed"
