#!/bin/bash
#SBATCH --job-name=catalog_build
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00

# catalog_build.sh — build a founder SNP catalog (PROPOSED parallel REFALT path).
#
# Part of an ALTERNATIVE ref/alt pipeline under evaluation. It does not touch the
# validated bam2bcf2REFALT.sh path. The idea: define the usable SNP set ONCE from
# the founders (a biological "every founder near-fixed and segregating" filter),
# instead of the current QUAL>59 filter that shifts with interval spec / cohort /
# BAQ. Samples are then counted against this fixed catalog (see catalog_count.sh).
#
# Steps:
#   1. Call variants from the founder BAMs only (BAQ off for determinism).
#   2. Drop SNPs within --snpgap bp of a founder indel (bcftools --SnpGap).
#   3. Keep biallelic SNPs where every founder has depth >= --min-dp and is
#      near-fixed (ALT freq <= --maxaf or >= 1-maxaf), and the site segregates
#      among founders (at least one REF-fixed and one ALT-fixed founder).
#   4. Emit catalog.tsv.gz (CHROM POS REF ALT), tabix-indexed, for mpileup -T.
#
# The fixation thresholds match the current pipeline's downstream good_SNPs step
# (REFALT2haps.code.R): near-fixed = ALT freq <= 0.03 or >= 0.97. The depth floor
# is stricter and explicit: every founder must have >= --min-dp (default 20x)
# reads, vs good_SNPs's any-depth. SNPs are catalogued genome-wide —
# heterochromatin is NOT censored ("you never know").
#
# TWO WAYS to specify founders:
#   (a) project    : --parfile hap_params.R --bamlist bam_list.txt
#                    Founder NAMES come from the config's `founders` vector
#                    (population founders + any tester strain); each is resolved
#                    to a BAM in bam_list by SM read-group tag. Use this per
#                    project (handles tester-strain designs).
#   (b) shared/direct: --founders <founder_bams.txt>
#                    Every BAM in the file is a founder (no name resolution). Use
#                    this to build the shared A-pop / B-pop catalogs once, e.g.
#                    from pipeline/helpfiles/{A,B}_founders.bams.txt.
#
# Usage:
#   sbatch catalog_build.sh --parfile helpfiles/<project>/hap_params.R \
#                           --bamlist helpfiles/<project>/bam_list.txt \
#                           --dir     process/<project>_catalog
#   sbatch catalog_build.sh --founders pipeline/helpfiles/B_founders.bams.txt \
#                           --dir     process/catalog_Bpop
#
# Options (defaults): --min-dp 20  --maxaf 0.03  --snpgap 5 (0 disables)
#                     --ref pipeline/ref/dm6.fa

set -euo pipefail

REF=pipeline/ref/dm6.fa
MIN_DP=20
MAXAF=0.03
SNPGAP=5
THREADS=8

while [[ $# -gt 0 ]]; do
  case $1 in
    --parfile)      PARFILE="$2";       shift 2 ;;
    --bamlist)      BAMLIST="$2";       shift 2 ;;
    --founders)     FOUNDER_BAMS="$2";  shift 2 ;;
    --dir)          DIR="$2";           shift 2 ;;
    --ref)          REF="$2";           shift 2 ;;
    --min-dp)       MIN_DP="$2";        shift 2 ;;
    --maxaf)        MAXAF="$2";         shift 2 ;;
    --snpgap)       SNPGAP="$2";        shift 2 ;;
    --threads)      THREADS="$2";       shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "${DIR:-}" ]] && { echo "Error: --dir required" >&2; exit 1; }

module load bcftools/1.21
module load samtools/1.10

mkdir -p "${DIR}"
raw="${DIR}/founders.calls.bcf"
cat="${DIR}/catalog.tsv.gz"
flist="${DIR}/founders.bams.txt"

# The catalog is built ONCE per project and reused. If it already exists (e.g.
# on a rerun to add samples), do not recall the founders. FORCE=1 rebuilds.
if [[ -s "$cat" && "${FORCE:-0}" != "1" ]]; then
  echo "catalog already exists ($cat); reusing (founders not recalled). FORCE=1 to rebuild."
  exit 0
fi

# Assemble the founder BAM list, either directly (shared) or by config (project).
if [[ -n "${FOUNDER_BAMS:-}" ]]; then
  grep -ve '^[[:space:]]*$' "$FOUNDER_BAMS" > "$flist"
  echo "founders (direct): $(grep -cve '^[[:space:]]*$' "$flist") BAM(s) from $FOUNDER_BAMS"
else
  [[ -z "${PARFILE:-}" ]] && { echo "Error: --parfile (with --bamlist) or --founders required" >&2; exit 1; }
  [[ -z "${BAMLIST:-}" ]] && { echo "Error: --bamlist required with --parfile" >&2; exit 1; }
  # Founder names, straight from the project config (hap_params.R `founders`).
  founders=$(grep -E '^[[:space:]]*founders' "$PARFILE" | grep -oE '"[^"]+"' | tr -d '"')
  [[ -z "$founders" ]] && { echo "Error: no 'founders' vector found in $PARFILE" >&2; exit 1; }
  # Resolve each founder name to a BAM in the bam list by its SM read-group tag.
  : > "$flist"
  while IFS= read -r bam; do
    [[ -z "$bam" ]] && continue
    sm=$(samtools view -H "$bam" | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print substr($i,4)}' | sort -u | head -1)
    for f in $founders; do
      if [[ "$sm" == "$f" ]]; then echo "$bam" >> "$flist"; break; fi
    done
  done < "$BAMLIST"
  nfound=$(grep -cve '^[[:space:]]*$' "$flist" || true)
  nwant=$(echo "$founders" | grep -cve '^[[:space:]]*$')
  echo "founders from $PARFILE ($nwant): $(echo $founders | tr '\n' ' ')"
  echo "matched $nfound founder BAM(s) by SM tag in $BAMLIST -> $flist"
  [[ "$nfound" -eq "$nwant" ]] || { echo "Error: matched $nfound of $nwant founders; check SM tags / bam_list" >&2; exit 1; }
fi

# 1-2. Call founders (BAQ off, -B), split multiallelics, mark SNPs near indels.
bcftools mpileup -B -q 20 -Q 20 --max-depth 1000 -a FORMAT/AD,FORMAT/DP \
    -f "$REF" --threads "$THREADS" -b "$flist" -Ou \
  | bcftools call -mv --threads "$THREADS" -Ou \
  | bcftools norm -f "$REF" -m - --threads "$THREADS" -Ou \
  | bcftools filter -g "$SNPGAP" --threads "$THREADS" -Ob > "$raw"

# 3. Founder-fixation filter on biallelic SNPs not flagged near an indel.
#    Genome-wide; heterochromatin is not censored.
#    -> catalog.bed: CHROM POS0 POS REF ALT   (POS is 1-based catalog position)
bcftools view -m2 -M2 -v snps -e 'FILTER~"SnpGap"' "$raw" \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
  | awk -F'\t' -v OFS='\t' -v mindp="$MIN_DP" -v maxaf="$MAXAF" '
    {
        segR = 0; segA = 0; ok = 1
        for (i = 5; i <= NF; i++) {            # one AD field per founder
            split($i, a, ",")
            r = a[1] + 0; v = a[2] + 0; dp = r + v
            if (dp < mindp) { ok = 0; break }  # every founder must be covered
            af = v / dp
            if (af > maxaf && af < 1 - maxaf) { ok = 0; break }   # not near-fixed
            if (af <= maxaf)     segR = 1
            if (af >= 1 - maxaf) segA = 1
        }
        if (ok && segR && segA) print $1, $2 - 1, $2, $3, $4     # segregating
    }' > "${DIR}/catalog.bed"

# 4. Positions file for `bcftools mpileup -T` (CHROM POS REF ALT), bgzipped+tabixed.
awk -v OFS='\t' '{print $1, $3, $4, $5}' "${DIR}/catalog.bed" \
  | sort -k1,1 -k2,2n | bgzip > "$cat"
tabix -s1 -b2 -e2 "$cat"

echo "catalog sites (genome-wide): $(wc -l < "${DIR}/catalog.bed")  -> $cat"
