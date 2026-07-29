#!/bin/bash
#SBATCH --job-name=catalog_filter
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:30:00

# catalog_filter.sh — apply thresholds to the annotated catalog -> catalog.tsv.gz.
#
# Part of the PROPOSED founder-catalog caller. Turns <catdir>/catalog.annot.tsv.gz
# (all candidate SNPs + annotations, from build_catalog.sh) into the bare
# CHROM POS REF,ALT positions file <catdir>/catalog.tsv.gz that call_samples.sh
# feeds to `bcftools -T`. Because the founder calling is NOT redone, changing a
# threshold (e.g. --snpgap 25 vs 5) is a fast re-cut, not an hours-long rebuild
# (XQTL2 #21). Also writes <catdir>/catalog.stats.txt (per-rule tally).
#
# A SNP is kept if, across the founders (an --exempt founder on the SNP's
# chromosome is skipped, as if not a founder):
#   snpgap    : >= --snpgap bp from the nearest founder indel
#   depth     : every founder has >= --min-dp reads
#   near-fixed: every founder's ALT freq <= --maxaf or >= 1-maxaf
#   segregating: >=1 founder REF-fixed and >=1 ALT-fixed
#
# Usage (standalone re-cut, or via build_catalog.sh):
#   sbatch catalog_filter.sh --catdir <catalog dir> \
#          [--min-dp 10 --maxaf 0.03 --snpgap 25 --exempt-founders B5:chr2L]

set -euo pipefail

MIN_DP=10
MAXAF=0.03
SNPGAP=20          # default from the AGE_SY indel-distance sweep (#28); tune per data
EXEMPT="B5:chr2L"

while [[ $# -gt 0 ]]; do
  case $1 in
    --catdir)          CATDIR="$2";  shift 2 ;;
    --min-dp)          MIN_DP="$2";  shift 2 ;;
    --maxaf)           MAXAF="$2";   shift 2 ;;
    --snpgap)          SNPGAP="$2";  shift 2 ;;
    --exempt-founders) EXEMPT="$2";  shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -z "${CATDIR:-}" ]] && { echo "Error: --catdir required" >&2; exit 1; }
annot="${CATDIR}/catalog.annot.tsv.gz"
[[ -f "$annot" ]] || { echo "Error: no catalog.annot.tsv.gz in $CATDIR (run build_catalog.sh first)" >&2; exit 1; }

module load bcftools/1.21
out="${CATDIR}/catalog.tsv.gz"
stats="${CATDIR}/catalog.stats.txt"
tmp="${CATDIR}/.catalog.stats.raw"

# Filter the annotated table. Kept SNPs -> stdout as CHROM POS REF,ALT; tally -> tmp.
# One POSITION at a time: the annotated table is coordinate-sorted, so a multiallelic
# founder site (>1 ALT, split by norm -m - into several rows at one POS) shows up as
# adjacent duplicate (CHROM,POS). Such positions are dropped (a pooled REF-vs-ALT
# counter is undefined there; XQTL2 #24). Single-biallelic positions go through the
# snpgap / depth / near-fixed / segregating filters.
zcat "$annot" \
  | awk -F'\t' -v OFS='\t' -v mindp="$MIN_DP" -v maxaf="$MAXAF" -v snpgap="$SNPGAP" \
        -v EXEMPT="$EXEMPT" -v statsfile="$tmp" '
    function decide(line,   nf, F, i, chr, segR, segA, nm, a, r, v, dp, af) {
        total++
        nf = split(line, F, "\t")
        chr = F[1]
        if ((F[5] + 0) < snpgap) { c_snpgap++; return }
        segR = 0; segA = 0
        for (i = 6; i <= nf; i++) {
            nm = fname[i]
            if (nm in exemptall)           continue
            if ((nm SUBSEP chr) in exempt) continue
            split(F[i], a, ","); r = a[1] + 0; v = a[2] + 0; dp = r + v
            if (dp < mindp)                   { c_depth++;    return }
            af = v / dp
            if (af > maxaf && af < 1 - maxaf) { c_notfixed++; return }
            if (af <= maxaf)     segR = 1
            if (af >= 1 - maxaf) segA = 1
        }
        if (!(segR && segA)) { c_notseg++; return }
        c_kept++
        print F[1], F[2], F[3]","F[4]
    }
    NR == 1 {
        for (i = 6; i <= NF; i++) { nm = $i; sub(/^AD_/, "", nm); fname[i] = nm }
        ne = split(EXEMPT, earr, ",")
        for (k = 1; k <= ne; k++) {
            e = earr[k]; if (e == "") continue
            ci = index(e, ":")
            if (ci > 0) exempt[substr(e,1,ci-1) SUBSEP substr(e,ci+1)] = 1
            else        exemptall[e] = 1
        }
        next
    }
    {
        key = $1 SUBSEP $2
        if (key != pk) {
            if (pk != "") { if (runcnt == 1) decide(buf); else { total++; c_multi++ } }
            pk = key; buf = $0; runcnt = 1
        } else runcnt++
    }
    END {
        if (pk != "") { if (runcnt == 1) decide(buf); else { total++; c_multi++ } }
        printf "candidates=%d multi_fail=%d snpgap_fail=%d depth_fail=%d notfixed_fail=%d notseg_fail=%d kept=%d\n",
               total, c_multi, c_snpgap, c_depth, c_notfixed, c_notseg, c_kept > statsfile
    }' \
  | sort -k1,1 -k2,2n | bgzip > "$out"
tabix -s1 -b2 -e2 "$out"

# Pretty-print the tally with the thresholds used.
awk -v mindp="$MIN_DP" -v maxaf="$MAXAF" -v snpgap="$SNPGAP" -v exempt="$EXEMPT" '
  { for (i=1;i<=NF;i++){ split($i,kv,"="); t[kv[1]]=kv[2] } }
  END {
    printf "catalog filter: min-dp=%s maxaf=%s snpgap=%s exempt=%s\n", mindp, maxaf, snpgap, exempt
    printf "  candidate positions      : %d\n", t["candidates"]
    printf "  dropped, multiallelic    : %d\n", t["multi_fail"]
    printf "  dropped, near indel      : %d\n", t["snpgap_fail"]
    printf "  dropped, founder depth   : %d\n", t["depth_fail"]
    printf "  dropped, not near-fixed  : %d\n", t["notfixed_fail"]
    printf "  dropped, not segregating : %d\n", t["notseg_fail"]
    printf "  KEPT (catalog)           : %d\n", t["kept"]
  }' "$tmp" > "$stats"
rm -f "$tmp"

echo "catalog: $(zcat "$out" | wc -l) sites -> $out"
cat "$stats"
