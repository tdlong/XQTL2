# catalog_merge.R — merge per-sample catalog counts into drop-in RefAlt.<chr>.txt
#
# Part of the PROPOSED parallel REFALT path (see catalog_build.sh / catalog_count.sh).
# Reads every counts/<sample>.tsv.gz, outer-joins them on the catalog site
# (CHROM,POS,REF,ALT), fills missing coverage with 0, and writes one
# RefAlt.<chr>.txt per chromosome in the SAME format the validated pipeline
# produces (CHROM POS REF_<name> ALT_<name> ...), so REFALT2haps runs unchanged.
#
# Column order is cosmetic here: REFALT2haps matches samples by column NAME, not
# position. REF_<name>/ALT_<name> are kept adjacent per sample for readability.
#
# Usage:
#   Rscript catalog_merge.R <output_dir> [chr1,chr2,...]

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args   <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
chrs   <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

files <- list.files(file.path(outdir, "counts"), pattern = "\\.tsv\\.gz$", full.names = TRUE)
if (length(files) == 0) stop("no counts/*.tsv.gz in ", outdir)

key <- c("CHROM", "POS", "REF", "ALT")
# read via gzip -dc so we don't depend on the R.utils package for .gz
tabs <- lapply(files, function(f) fread(cmd = paste("gzip -dc", shQuote(f)), sep = "\t", header = TRUE))
merged <- reduce(tabs, full_join, by = key)

# missing (sample, site) coverage -> 0 reads
countcols <- setdiff(names(merged), key)
setDT(merged)
merged[, (countcols) := lapply(.SD, function(x) fifelse(is.na(x), 0L, as.integer(x))), .SDcols = countcols]

# drop the allele columns; RefAlt.<chr>.txt is CHROM POS + REF_/ALT_ per sample
merged[, c("REF", "ALT") := NULL]
setorder(merged, CHROM, POS)

for (chr in chrs) {
  sub <- merged[CHROM == chr]
  f <- file.path(outdir, paste0("RefAlt.", chr, ".txt"))
  fwrite(sub, f, sep = "\t")
  cat(sprintf("wrote %s (%d sites, %d samples)\n", f, nrow(sub), length(countcols) / 2))
}
