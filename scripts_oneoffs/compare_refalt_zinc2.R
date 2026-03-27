#!/usr/bin/env Rscript
# Validate the founder-sites RefAlt approach against the original joint-called RefAlt.
# Comparison is made AFTER applying the good_SNPs filter from REFALT2haps.code.R,
# since that filter — not raw bcftools output — defines the sites actually used downstream.
#
# Usage: Rscript compare_refalt_zinc2.R <original_dir> <new_dir>
#   original_dir : process/ZINC2/          (joint-called RefAlt)
#   new_dir      : process/ZINC2_fromsites/ (founder-sites RefAlt)

library(data.table)
library(tidyverse)

args         <- commandArgs(trailingOnly = TRUE)
original_dir <- args[1]
new_dir      <- args[2]

chrs <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

# B-population founder SM tag names (as they appear in RefAlt column headers)
founders <- c("AB8", "B1", "B2", "B3", "B4", "B5", "B6", "B7")

# ---------------------------------------------------------------------------
# Apply the good_SNPs filter from REFALT2haps.code.R to a wide RefAlt table.
# Returns a data.table of CHROM+POS for sites passing all three criteria:
#   zeros==0      : every founder has reads (N > 0)
#   notfixed==0   : no founder appears heterozygous (freq between 0.03-0.97)
#   informative   : site is polymorphic between founders
# ---------------------------------------------------------------------------
good_snps_filter <- function(dt, founders) {
  # Extract founder REF/ALT columns and pivot to long
  founder_ref <- paste0("REF_", founders)
  founder_alt <- paste0("ALT_", founders)
  present <- founders[founder_ref %in% names(dt) & founder_alt %in% names(dt)]

  if (length(present) == 0) stop("No founder columns found in RefAlt table")

  # Build long founder table: one row per site per founder
  long <- rbindlist(lapply(present, function(f) {
    dt[, .(CHROM, POS,
           founder = f,
           REF     = get(paste0("REF_", f)),
           ALT     = get(paste0("ALT_", f)))]
  }))
  long[, N    := REF + ALT]
  long[, freq := ifelse(N == 0, NA_real_, REF / N)]

  long[, .(
    zeros       = sum(N == 0),
    notfixed    = sum(N != 0 & freq > 0.03 & freq < 0.97),
    informative = (sum(freq, na.rm=TRUE) > 0.05 | sum(freq, na.rm=TRUE) < (length(present) * 0.95))
  ), by = .(CHROM, POS)][zeros == 0 & notfixed == 0 & informative == TRUE, .(CHROM, POS)]
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
results <- lapply(chrs, function(chr) {

  orig_file <- file.path(original_dir, paste0("RefAlt.", chr, ".txt"))
  new_file  <- file.path(new_dir,      paste0("RefAlt.", chr, ".txt"))

  if (!file.exists(orig_file)) { message("MISSING original: ", orig_file); return(NULL) }
  if (!file.exists(new_file))  { message("MISSING new: ",      new_file);  return(NULL) }

  message("\n=== ", chr, " ===")

  # Original uses tabs in header, spaces in data — read.table handles mixed whitespace
  orig <- as.data.table(read.table(orig_file, header=TRUE, check.names=FALSE))
  new  <- fread(new_file)
  setkey(orig, CHROM, POS)
  setkey(new,  CHROM, POS)

  message(sprintf("  Raw sites — original: %d   new: %d", nrow(orig), nrow(new)))

  # --- Apply good_SNPs filter to both ---
  good_orig <- good_snps_filter(orig, founders)
  good_new  <- good_snps_filter(new,  founders)
  setkey(good_orig, CHROM, POS)
  setkey(good_new,  CHROM, POS)

  n_good_orig <- nrow(good_orig)
  n_good_new  <- nrow(good_new)
  shared      <- fintersect(good_orig, good_new)
  only_orig   <- fsetdiff(good_orig, good_new)
  only_new    <- fsetdiff(good_new,  good_orig)

  message(sprintf("  good_SNPs sites — original: %d   new: %d", n_good_orig, n_good_new))
  message(sprintf("  Shared good_SNPs: %d  (%.2f%% of original)", nrow(shared), 100*nrow(shared)/n_good_orig))
  message(sprintf("  Only in original: %d   Only in new: %d", nrow(only_orig), nrow(only_new)))

  if (nrow(only_new) > 0) {
    message("  *** WARNING: good_SNPs in new but not original — check founder catalog ***")
    print(only_new[1:min(5, nrow(only_new))])
  }

  # --- Exact match of pool counts at shared good_SNPs sites ---
  get_pool_samples <- function(dt) {
    cols <- setdiff(names(dt), c("CHROM","POS"))
    cols <- cols[!grepl(paste(paste0("_(", paste(founders, collapse="|"), ")$"), sep=""), cols)]
    unique(sub("^(REF|ALT)_", "", cols))
  }

  common_pools <- intersect(get_pool_samples(orig), get_pool_samples(new))
  message(sprintf("  Common pool samples: %d", length(common_pools)))

  orig_at_shared <- orig[shared, on = c("CHROM","POS")]
  new_at_shared  <- new[shared,  on = c("CHROM","POS")]

  mismatches <- rbindlist(lapply(common_pools, function(s) {
    rc <- paste0("REF_", s); ac <- paste0("ALT_", s)
    if (!rc %in% names(orig_at_shared) || !rc %in% names(new_at_shared)) return(NULL)
    n_ref <- sum(orig_at_shared[[rc]] != new_at_shared[[rc]])
    n_alt <- sum(orig_at_shared[[ac]] != new_at_shared[[ac]])
    if (n_ref > 0 || n_alt > 0)
      data.table(sample=s, n_ref_diff=n_ref, n_alt_diff=n_alt)
  }))

  if (nrow(mismatches) == 0) {
    message(sprintf("  EXACT MATCH: all %d pool samples identical at %d shared good_SNPs sites",
                    length(common_pools), nrow(shared)))
  } else {
    message(sprintf("  *** COUNT MISMATCH in %d samples at shared good_SNPs sites ***", nrow(mismatches)))
    print(mismatches)
  }

  list(chr=chr,
       raw_orig=nrow(orig), raw_new=nrow(new),
       good_orig=n_good_orig, good_new=n_good_new,
       n_shared=nrow(shared), n_only_orig=nrow(only_orig), n_only_new=nrow(only_new),
       n_mismatched_samples=nrow(mismatches))
})

# --- Summary ---
message("\n\n=== SUMMARY ===")
summary_dt <- rbindlist(lapply(results, function(r) {
  if (is.null(r)) return(NULL)
  data.table(chr          = r$chr,
             raw_orig     = r$raw_orig,
             raw_new      = r$raw_new,
             good_orig    = r$good_orig,
             good_new     = r$good_new,
             shared       = r$n_shared,
             pct_shared   = round(100 * r$n_shared / r$good_orig, 2),
             only_orig    = r$n_only_orig,
             only_new     = r$n_only_new,
             pool_mismatch= r$n_mismatched_samples)
}))
print(summary_dt)
