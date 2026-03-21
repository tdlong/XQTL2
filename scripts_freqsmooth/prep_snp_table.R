#!/usr/bin/env Rscript
# prep_snp_table.R
#
# Extract founder columns from FREQ_SNPs.cM.txt.gz, discarding the pool
# frequency columns (~28 cols) not used by snp_scan.R.  Writes two files:
#
#   FREQ_SNPs_Apop.cM.txt.gz  —  CHROM POS A1 A2 A3 A4 A5 A6 A7 AB8       cM
#   FREQ_SNPs_Bpop.cM.txt.gz  —  CHROM POS                  AB8 B1..B7    cM
#
# AB8 is shared between A and B populations and is included in both files.
# Column order within each file matches the source file (A1..AB8, AB8..B7).
#
# Usage:
#   Rscript scripts_freqsmooth/prep_snp_table.R /path/to/FREQ_SNPs.cM.txt.gz
#
# Output files are written to the same directory as the input.

suppressPackageStartupMessages(library(data.table))

args   <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
if (is.na(infile)) stop("Usage: Rscript prep_snp_table.R <FREQ_SNPs.cM.txt.gz>")
outdir <- dirname(infile)

apop_cols <- c("CHROM","POS","A1","A2","A3","A4","A5","A6","A7","AB8","cM")
bpop_cols <- c("CHROM","POS","AB8","B1","B2","B3","B4","B5","B6","B7","cM")
keep      <- union(apop_cols, bpop_cols)

cat("Reading", infile, "...\n")
snp <- fread(infile, select = keep)
cat(sprintf("  %d SNPs read\n", nrow(snp)))

out_a <- file.path(outdir, "FREQ_SNPs_Apop.cM.txt.gz")
out_b <- file.path(outdir, "FREQ_SNPs_Bpop.cM.txt.gz")

cat("Writing", out_a, "...\n")
fwrite(snp[, ..apop_cols], out_a, sep = "\t", compress = "gzip")

cat("Writing", out_b, "...\n")
fwrite(snp[, ..bpop_cols], out_b, sep = "\t", compress = "gzip")

cat("Done.\n")
