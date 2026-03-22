#!/bin/bash
# concat_snp_scans.sh  <outdir>
#
# Merges per-chromosome snp_scan files into a single <name>.snp_scan.txt.
# Usage:
#   bash scripts_freqsmooth/concat_snp_scans.sh process/proj/SCAN_NAME

module load R/4.2.2

mydir=$1
name=$(basename "$mydir")

Rscript - <<REOF
data_path <- "$mydir"
name      <- "$name"
library(tidyverse)
files <- dir(data_path, pattern = paste0(name, "\\\\.snp_scan\\\\..*\\\\.txt"), full.names = TRUE)
files <- grep("chr", files, value = TRUE)
cat("Merging", length(files), "SNP scan files\\n")
df <- files %>%
  map(~ as_tibble(read.table(.x, header = TRUE))) %>%
  bind_rows()
outfile <- file.path(data_path, paste0(name, ".snp_scan.txt"))
write.table(df, outfile, quote = FALSE)
cat("Written:", outfile, "\\n")
cat(sprintf("Total SNPs: %d\\n", nrow(df)))
REOF
