# catalog_merge.R — merge per-sample catalog counts into drop-in RefAlt.<chr>.txt
#
# Part of the PROPOSED founder-catalog caller (a worker driven by call_samples.sh).
# Reads every <calls dir>/counts/<sample>.tsv.gz (columns CHROM POS REF_<name>
# ALT_<name>), joins them on (CHROM,POS) — the catalog already fixes the alleles,
# so there are no allele columns to key on — fills missing coverage with 0, and
# writes one <calls dir>/RefAlt.<chr>.txt per chromosome, in the SAME format the
# validated pipeline produces (CHROM POS REF_<name> ALT_<name> ...) so REFALT2haps
# runs unchanged.
#
# Memory: all count files are the same catalog sites in the same order, so this
# just lines up REF/ALT columns. <2M SNPs x ~100 samples x 8 bytes ~= 1.6GB;
# samples are joined one at a time (peak ~one copy, ~3GB). Standard 6G/core ample.
#
# Usage:
#   Rscript catalog_merge.R <calls dir> [chr1,chr2,...]

suppressPackageStartupMessages(library(data.table))

args     <- commandArgs(trailingOnly = TRUE)
callsdir <- args[1]
chrs     <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

files <- list.files(file.path(callsdir, "counts"), pattern = "\\.tsv\\.gz$", full.names = TRUE)
if (length(files) == 0) stop("no counts/*.tsv.gz in ", callsdir)
key <- c("CHROM", "POS")

# Join one sample at a time; peak memory stays ~one copy of the growing table.
merged <- NULL
for (f in files) {
  dt <- fread(cmd = paste("gzip -dc", shQuote(f)), sep = "\t", header = TRUE)   # CHROM POS REF_<name> ALT_<name>
  merged <- if (is.null(merged)) dt else merge(merged, dt, by = key, all = TRUE)
}

cc <- setdiff(names(merged), key)
merged[, (cc) := lapply(.SD, function(x) fifelse(is.na(x), 0L, as.integer(x))), .SDcols = cc]
setorder(merged, CHROM, POS)

for (chr in chrs) {
  sub <- merged[CHROM == chr]
  f <- file.path(callsdir, paste0("RefAlt.", chr, ".txt"))
  fwrite(sub, f, sep = "\t")
  cat(sprintf("wrote %s (%d sites, %d samples)\n", f, nrow(sub), length(cc) / 2))
}
