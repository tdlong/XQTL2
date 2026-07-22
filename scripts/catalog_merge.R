# catalog_merge.R — merge per-sample catalog counts into drop-in RefAlt.<chr>.txt
#
# Part of the PROPOSED parallel REFALT path (see catalog_count.sh). Reads every
# counts/<sample>.tsv.gz, joins them on the catalog site (CHROM,POS,REF,ALT),
# fills missing coverage with 0, and writes one RefAlt.<chr>.txt per chromosome
# in the SAME format the validated pipeline produces (CHROM POS REF_<name>
# ALT_<name> ...), so REFALT2haps runs unchanged.
#
# Memory: all count files are the SAME catalog sites in the SAME order (each is
# mpileup -T against the one catalog), so this is really just lining up REF/ALT
# columns. The whole-genome catalog is <2M SNPs; each sample adds two 4-byte
# integers = 8 bytes/site (~16 MB/sample), so the full wide table for ~100
# samples is ~1.6GB. Samples are joined one at a time, keeping peak to ~one copy
# (~3GB). Standard 6G/core is ample; no highmem, no chromosome subsetting.
#
# Usage:
#   Rscript catalog_merge.R <output_dir> [chr1,chr2,...]

suppressPackageStartupMessages(library(data.table))

args   <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
chrs   <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

files <- list.files(file.path(outdir, "counts"), pattern = "\\.tsv\\.gz$", full.names = TRUE)
if (length(files) == 0) stop("no counts/*.tsv.gz in ", outdir)
key <- c("CHROM", "POS", "REF", "ALT")

# Join one sample at a time; peak memory stays ~one copy of the growing table.
merged <- NULL
for (f in files) {
  dt <- fread(cmd = paste("gzip -dc", shQuote(f)), sep = "\t", header = TRUE)   # CHROM POS REF ALT REF_<name> ALT_<name>
  merged <- if (is.null(merged)) dt else merge(merged, dt, by = key, all = TRUE)
}

cc <- setdiff(names(merged), key)
merged[, (cc) := lapply(.SD, function(x) fifelse(is.na(x), 0L, as.integer(x))), .SDcols = cc]
merged[, c("REF", "ALT") := NULL]      # RefAlt.<chr>.txt is CHROM POS + REF_/ALT_ per sample
setorder(merged, CHROM, POS)

for (chr in chrs) {
  sub <- merged[CHROM == chr]
  f <- file.path(outdir, paste0("RefAlt.", chr, ".txt"))
  fwrite(sub, f, sep = "\t")
  cat(sprintf("wrote %s (%d sites, %d samples)\n", f, nrow(sub), length(cc) / 2))
}
