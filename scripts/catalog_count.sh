#!/bin/bash
#SBATCH --job-name=catalog_count
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-00:00:00
# Submitted by run_refalt.catalog.sh as an array with ONE SAMPLE per task
# (--array=1-<#BAMs>). Counting a fixed catalog is not split by chromosome —
# one whole-genome job per sample. Since that is a pile of SLURM jobs, `seff`
# them afterward for real per-sample timings.

# catalog_count.sh — count REF/ALT reads for ONE sample against the catalog.
#
# The output directory is PERSISTENT project state: a sample already counted
# (its counts/<name>.tsv.gz exists) is SKIPPED, so re-running after appending new
# samples to the bam list only counts the new ones — the founder catalog is never
# rebuilt and prior samples are never recounted. Counting is deterministic (BAQ
# off -B; positions fixed by -T catalog), so counts do not depend on interval
# spec or on which other samples are in the run.
#
# Column name comes from the BAM's SM read-group tag, matching the validated
# pipeline's naming so the merge is drop-in. Set FORCE=1 to recount.
#
# Usage (via run_refalt.catalog.sh):
#   sbatch --array=1-<#BAMs> catalog_count.sh <bam_list> <output_dir>

set -euo pipefail

ref="pipeline/ref/dm6.fa"
bams=$1
output=$2
cat="${output}/catalog.tsv.gz"

[[ -f "$cat" ]] || { echo "Error: catalog not found: $cat (run catalog_build.sh + catalog_gather.sh first)" >&2; exit 1; }

# NB: bcftools/1.21 and samtools/1.10 link incompatible htslib versions and must
# never be loaded in the same shell (bcftools then fails with exit 127). Load
# bcftools here; read the SM tag with samtools in an isolated subshell below.
module load bcftools/1.21
mkdir -p "${output}/counts"

bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$bams")
[[ -z "$bam" ]] && { echo "No BAM at line $SLURM_ARRAY_TASK_ID of $bams" >&2; exit 1; }

# Resolve the sample name (SM tag) up front so an already-counted sample is
# skipped WITHOUT redoing the pileup — this is what makes adding samples cheap.
# samtools runs in an isolated subshell so its htslib never enters this shell's
# bcftools environment.
name=$( module load samtools/1.10; samtools view -H "$bam" | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print substr($i,4)}' | sort -u | head -1 )
[[ -z "$name" ]] && { echo "No SM tag in $bam" >&2; exit 1; }
out="${output}/counts/${name}.tsv.gz"

if [[ -s "$out" && "${FORCE:-0}" != "1" ]]; then
  echo "skip ${name}: already counted (${out}); set FORCE=1 to recount"
  exit 0
fi

# Count reads for the CATALOG's fixed REF/ALT at each site. mpileup reports ALT
# per sample (<*>, C,<*>, ...), which would make AD{1} mean different alleles in
# different samples and fragment sites at merge time (issue #14). `call -m -C
# alleles -T catalog` constrains genotyping to the catalog's REF/ALT, so AD{0}=
# REF reads and AD{1}=ALT reads mean the same allele in every sample. -C alleles
# also emits every catalog site (0 coverage -> AD 0,0), so the tables align.
tmp="${output}/counts/tmp.${SLURM_ARRAY_TASK_ID}.bcf"
bcftools mpileup -B -q 20 -Q 20 --max-depth 2000 -T "$cat" -a FORMAT/AD \
    -f "$ref" "$bam" -Ou \
  | bcftools call -m -C alleles -T "$cat" -Ob > "$tmp"

# Emit CHROM POS REF_<name> ALT_<name>, keyed on position (alleles fixed by catalog).
# Not tabix-indexed: catalog_merge.R reads each file whole via `gzip -dc`, and the
# header line ("CHROM POS ...") is not a genomic interval, so indexing it is both
# pointless and an error.
{
  printf 'CHROM\tPOS\tREF_%s\tALT_%s\n' "$name" "$name"
  bcftools query -f '%CHROM\t%POS[\t%AD]\n' "$tmp" \
    | awk -F'\t' -v OFS='\t' '{ split($3,a,","); print $1,$2,(a[1]==""||a[1]=="."?0:a[1]),(a[2]==""||a[2]=="."?0:a[2]) }'
} | bgzip > "$out"

rm -f "$tmp"
echo "counted ${name} -> ${out}"
