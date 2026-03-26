#!/usr/bin/env Rscript
# Consolidate per-sample RefAlt files into a single wide RefAlt table.
#
# Usage: Rscript consolidate_refalt_fromsites.R <sample_dir> <chr> <sites_vcf> <outfile>
#   sample_dir : directory containing RefAlt.<samplename>.<chr>.txt files
#   chr        : chromosome (e.g. chrX)
#   sites_vcf  : founder_sites.<chr>.vcf.gz — defines the complete site list
#   outfile    : output path for wide RefAlt.<chr>.txt
#
# Sites present in the founder catalog but absent from a sample file (zero coverage
# or zero alt reads) are filled with 0 for both REF and ALT counts.
#
# To add new samples to an existing project:
#   1. Run bam2bcf2REFALT_fromsites.sh for each new sample
#   2. Re-run this script with all per-sample files (old + new) in sample_dir

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
sample_dir <- args[1]
chr        <- args[2]
sites_vcf  <- args[3]
outfile    <- args[4]

# --- Load founder site catalog as the master row set ---
# Read VCF skipping comment lines; we only need CHROM and POS
sites_raw <- fread(cmd = paste("bcftools query -f'%CHROM\\t%POS\\n'", sites_vcf),
                   col.names = c("CHROM", "POS"))
setkey(sites_raw, CHROM, POS)

# --- Load all per-sample RefAlt files for this chromosome ---
pattern <- paste0("RefAlt\\.", "*.", chr, "\\.txt$")
files   <- list.files(sample_dir, pattern = paste0(chr, "\\.txt$"), full.names = TRUE)
files   <- files[grepl("^RefAlt\\.", basename(files))]

if (length(files) == 0) stop("No per-sample RefAlt files found in: ", sample_dir)

message("Found ", length(files), " sample files for ", chr)

sample_tables <- lapply(files, function(f) {
  dt <- fread(f, sep = "\t")
  setkey(dt, CHROM, POS)
  dt
})

# --- Left-join each sample onto the master site list, fill missing with 0 ---
result <- copy(sites_raw)

for (dt in sample_tables) {
  # Get the REF_ and ALT_ column names for this sample
  sample_cols <- setdiff(names(dt), c("CHROM", "POS"))
  ref_col <- sample_cols[1]
  alt_col <- sample_cols[2]

  result <- merge(result, dt, by = c("CHROM", "POS"), all.x = TRUE)
  result[is.na(get(ref_col)), (ref_col) := 0L]
  result[is.na(get(alt_col)), (alt_col) := 0L]
}

setorder(result, CHROM, POS)

message("Writing ", nrow(result), " sites x ", ncol(result), " columns to ", outfile)
fwrite(result, outfile, sep = "\t")
