#!/bin/bash
#SBATCH --job-name=catalog_build
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00
# Submitted by run_refalt.catalog.sh as a 5-task array, ONE CHROMOSOME per task
# (bcftools mpileup is single-threaded, so per-chromosome parallelism matters —
# same reason bam2bcf2REFALT.sh is a per-chromosome array).

# catalog_build.sh — build a founder SNP catalog, one chromosome per array task.
#
# Part of the PROPOSED parallel REFALT path (not the validated bam2bcf2REFALT.sh).
# Defines the usable SNP set ONCE from the founders (every founder near-fixed and
# segregating), instead of the current QUAL>59 filter. Samples are then counted
# against this fixed catalog (catalog_count.sh).
#
# Per task (chromosome <chr>):
#   1. Call the founder BAMs on <chr> (BAQ off for determinism).
#   2. Drop SNPs within --snpgap bp of a founder indel (bcftools --SnpGap).
#   3. Keep biallelic SNPs where every founder has depth >= --min-dp and is
#      near-fixed (ALT freq <= --maxaf or >= 1-maxaf), and the site segregates.
#   4. Emit catalog.<chr>.bed for this chromosome. catalog_gather.sh then
#      concatenates the per-chromosome pieces into one genome-wide catalog.tsv.gz.
#
# Thresholds match the current good_SNPs step (near-fixed = freq <= 0.03 or
# >= 0.97); the depth floor is an explicit --min-dp (default 10x, required of
# every founder — some B-pop founders, e.g. the reference-subtracted B5, run
# shallow, so 20x across all founders throws out most sites).
# Catalogued genome-wide within each chromosome — heterochromatin not censored.
#
# TWO WAYS to specify founders:
#   (a) project     : --parfile hap_params.R --bamlist bam_list.txt
#                     Founder NAMES from the config's `founders` vector, each
#                     resolved to a BAM in bam_list by SM read-group tag.
#   (b) shared/direct: --founders <founder_bams.txt>  (all BAMs are founders)
#
# Usage (via run_refalt.catalog.sh, which sets --array=1-5):
#   sbatch --array=1-5 catalog_build.sh --parfile <hap_params.R> --bamlist <bam_list> --dir <dir>
#   sbatch --array=1-5 catalog_build.sh --founders <founder_bams.txt> --dir <dir>
#
# Options (defaults): --min-dp 10  --maxaf 0.03  --snpgap 5 (0 disables)
#                     --exempt-founders B5:chr2L  --ref pipeline/ref/dm6.fa
# --exempt-founders: comma-separated founders dropped from the rules (still counted
# as samples). Each is NAME (all chromosomes) or NAME:CHR (that chromosome only).
# Default B5:chr2L (B5's chr2L is a shallow, circular reconstructed haplotype; B5 is
# a normal founder elsewhere). Harmless if the name is absent (e.g. A-pop).

set -euo pipefail

REF=pipeline/ref/dm6.fa
MIN_DP=10
MAXAF=0.03
SNPGAP=5
THREADS=2   # bcftools mpileup pileup is single-threaded; --threads only helps BGZF I/O
EXEMPT="B5:chr2L" # founders dropped from the rule evaluation (as if not a founder) —
            # but still counted as samples in RefAlt. Comma-separated; each entry is
            # NAME (all chromosomes) or NAME:CHR (that chromosome only). Default
            # B5:chr2L: B5's chr2L is a reconstructed haplotype (see data/founders/
            # FOUNDERS.md) — shallow AND its alleles are circular there — but B5 is a
            # normal ~15-22x founder on the other arms, so it is only exempt on chr2L.

while [[ $# -gt 0 ]]; do
  case $1 in
    --parfile)          PARFILE="$2";       shift 2 ;;
    --bamlist)          BAMLIST="$2";       shift 2 ;;
    --founders)         FOUNDER_BAMS="$2";  shift 2 ;;
    --dir)              DIR="$2";           shift 2 ;;
    --ref)              REF="$2";           shift 2 ;;
    --min-dp)           MIN_DP="$2";        shift 2 ;;
    --maxaf)            MAXAF="$2";         shift 2 ;;
    --snpgap)           SNPGAP="$2";        shift 2 ;;
    --exempt-founders)  EXEMPT="$2";        shift 2 ;;
    --threads)          THREADS="$2";       shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "${DIR:-}" ]] && { echo "Error: --dir required" >&2; exit 1; }

chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$((SLURM_ARRAY_TASK_ID - 1))]:-}
[[ -z "$mychr" ]] && { echo "No chromosome for array index ${SLURM_ARRAY_TASK_ID:-<unset>}" >&2; exit 1; }

# NB: bcftools/1.21 and samtools/1.10 link incompatible htslib versions and must
# never be loaded in the same shell (bcftools then fails the htslib symbol lookup,
# exit 127). Load bcftools here; read SM tags with samtools in an isolated subshell.
module load bcftools/1.21

mkdir -p "${DIR}"
raw="${DIR}/founders.calls.${mychr}.bcf"
bed="${DIR}/catalog.${mychr}.bed"
flist="${DIR}/founders.bams.${mychr}.txt"   # per-chr, so array tasks never race

# The catalog is built ONCE per project and reused. If this chromosome's catalog
# piece already exists (e.g. on a rerun to add samples), do not recall the founders.
if [[ -s "$bed" && "${FORCE:-0}" != "1" ]]; then
  echo "${mychr}: catalog piece exists ($bed); reusing (founders not recalled). FORCE=1 to rebuild."
  exit 0
fi

# Assemble the founder BAM list, either directly (shared) or by config (project).
if [[ -n "${FOUNDER_BAMS:-}" ]]; then
  grep -ve '^[[:space:]]*$' "$FOUNDER_BAMS" > "$flist"
  echo "${mychr}: founders (direct): $(grep -cve '^[[:space:]]*$' "$flist") BAM(s)"
else
  [[ -z "${PARFILE:-}" ]] && { echo "Error: --parfile (with --bamlist) or --founders required" >&2; exit 1; }
  [[ -z "${BAMLIST:-}" ]] && { echo "Error: --bamlist required with --parfile" >&2; exit 1; }
  founders=$(grep -E '^[[:space:]]*founders' "$PARFILE" | grep -oE '"[^"]+"' | tr -d '"')
  [[ -z "$founders" ]] && { echo "Error: no 'founders' vector found in $PARFILE" >&2; exit 1; }
  # samtools loaded ONLY inside this subshell, so its htslib never enters the
  # bcftools environment of the parent shell (see the module-load note above).
  ( module load samtools/1.10
    while IFS= read -r bam; do
      [[ -z "$bam" ]] && continue
      sm=$(samtools view -H "$bam" | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print substr($i,4)}' | sort -u | head -1)
      for f in $founders; do
        if [[ "$sm" == "$f" ]]; then echo "$bam"; break; fi
      done
    done < "$BAMLIST"
  ) > "$flist"
  nfound=$(grep -cve '^[[:space:]]*$' "$flist" || true)
  nwant=$(echo "$founders" | grep -cve '^[[:space:]]*$')
  echo "${mychr}: matched $nfound of $nwant founder BAM(s) by SM tag"
  [[ "$nfound" -eq "$nwant" ]] || { echo "Error: matched $nfound of $nwant founders; check SM tags / bam_list" >&2; exit 1; }
fi

# 1-2. Call founders on this chromosome (BAQ off), split multiallelics, mark SNPs near indels.
bcftools mpileup -B -q 20 -Q 20 --max-depth 1000 -a FORMAT/AD,FORMAT/DP \
    -r "$mychr" -f "$REF" --threads "$THREADS" -b "$flist" -Ou \
  | bcftools call -mv --threads "$THREADS" -Ou \
  | bcftools norm -f "$REF" -m - --threads "$THREADS" -Ou \
  | bcftools filter -g "$SNPGAP" --threads "$THREADS" -Ob > "$raw"

# Founder order in the BCF == the AD column order below (query field 4+i for the
# i-th founder). Turn --exempt-founders names into the AD fields to skip.
fnames=( $(bcftools query -l "$raw") )
skipcols=""
for i in "${!fnames[@]}"; do
  for e in ${EXEMPT//,/ }; do
    ename="${e%%:*}"                       # "B5"    from "B5:chr2L"
    echr=""; [[ "$e" == *:* ]] && echr="${e#*:}"   # "chr2L" (empty = all chromosomes)
    if [[ "${fnames[$i]}" == "$ename" ]] && { [[ -z "$echr" ]] || [[ "$echr" == "$mychr" ]]; }; then
      skipcols="${skipcols} $((5 + i))"
    fi
  done
done
[[ -n "${skipcols// }" ]] && echo "${mychr}: exempting founder field(s)${skipcols} (${EXEMPT}) from depth/fixation gates"

# 3. Founder-fixation filter -> catalog.<chr>.bed: CHROM POS0 POS REF ALT
#    Exempt founders (e.g. shallow B5) are skipped in the loop: they neither gate
#    on depth nor contribute to the near-fixed / segregation tests.
bcftools view -m2 -M2 -v snps -e 'FILTER~"SnpGap"' "$raw" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
  | awk -F'\t' -v OFS='\t' -v mindp="$MIN_DP" -v maxaf="$MAXAF" -v skipcols="$skipcols" '
    BEGIN { n = split(skipcols, s, " "); for (k = 1; k <= n; k++) if (s[k] != "") skip[s[k]] = 1 }
    {
        segR = 0; segA = 0; ok = 1
        for (i = 5; i <= NF; i++) {
            if (i in skip) continue                 # exempt founder: free pass
            split($i, a, ",")
            r = a[1] + 0; v = a[2] + 0; dp = r + v
            if (dp < mindp) { ok = 0; break }
            af = v / dp
            if (af > maxaf && af < 1 - maxaf) { ok = 0; break }
            if (af <= maxaf)     segR = 1
            if (af >= 1 - maxaf) segA = 1
        }
        if (ok && segR && segA) print $1, $2 - 1, $2, $3, $4
    }' > "$bed"

echo "${mychr}: catalog sites $(wc -l < "$bed")  -> $bed (catalog_gather.sh assembles catalog.tsv.gz)"
