#!/bin/bash
# concat_scans.sh — merge per-chromosome scan files
#
# Usage:
#   bash scripts/concat_scans.sh <scan_dir>           # hap scan only
#   bash scripts/concat_scans.sh --snp <scan_dir>     # SNP scan only
#
# Hap scan mode (default):
#   Merges .scan.<chr>.txt and .meansBySample.<chr>.txt, generates Manhattan
#   plots, and bundles everything into a tarball.
#
# SNP scan mode (--snp):
#   Merges .snp_scan.<chr>.txt and .snp_meansBySample.<chr>.txt into single files.

module load R/4.2.2

MODE=hap
if [[ "$1" == "--snp" ]]; then
    MODE=snp
    shift
fi

mydir=$1
name=$(basename "$mydir")

if [[ "$MODE" == "hap" ]]; then
    Rscript scripts/concat_Chromosome_Scans.R "$mydir"
    tar -czvf "$mydir/${name}.tar.gz" -C "$mydir" \
        "${name}.scan.txt" \
        "${name}.meansBySample.txt" \
        "${name}.5panel.cM.png" \
        "${name}.5panel.Mb.png" \
        "${name}.Manhattan.png"
else
    # merge snp_scan (small — read into R)
    Rscript - <<REOF
data_path <- "$mydir"
name      <- "$name"
library(tidyverse)
files <- dir(data_path, pattern = paste0(name, "\\\\.snp_scan\\\\..*\\\\.txt"), full.names = TRUE)
files <- grep("chr", files, value = TRUE)
cat("Merging", length(files), "SNP scan files\n")
df <- files %>%
  map(~ as_tibble(read.table(.x, header = TRUE))) %>%
  bind_rows()
outfile <- file.path(data_path, paste0(name, ".snp_scan.txt"))
write.table(df, outfile, quote = FALSE)
cat("Written:", outfile, "\n")
cat(sprintf("Total SNPs: %d\n", nrow(df)))
REOF

    # merge snp_meansBySample (large — stream with head/tail, no R needed)
    mfiles=($(ls "$mydir"/${name}.snp_meansBySample.chr*.txt 2>/dev/null | sort))
    moutfile="$mydir/${name}.snp_meansBySample.txt"
    if [[ ${#mfiles[@]} -gt 0 ]]; then
        echo "Merging ${#mfiles[@]} SNP meansBySample files"
        head -1 "${mfiles[0]}" > "$moutfile"
        for f in "${mfiles[@]}"; do
            tail -n +2 "$f" >> "$moutfile"
        done
        echo "Written: $moutfile ($(wc -l < "$moutfile") lines)"
    fi
fi
