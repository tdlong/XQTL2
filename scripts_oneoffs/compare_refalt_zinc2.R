#!/usr/bin/env Rscript
# Validate the founder-sites RefAlt approach against the original joint-called RefAlt.
#
# Usage: Rscript compare_refalt_zinc2.R <original_dir> <new_dir>
#   original_dir : process/ZINC2/         (joint-called RefAlt files)
#   new_dir      : process/ZINC2_fromsites/ (consolidated founder-sites RefAlt files)
#
# What this checks:
#   1. Site coverage  — how many sites in each, how many shared, unique to each
#   2. Exact match    — at shared sites, counts must be bit-for-bit identical
#                       ANY mismatch is a red flag that needs investigation
#   3. Dropped sites  — sites in original but not new: are they near the QUAL threshold?
#                       Are they driven by pools (not founders)?
#   4. Gained sites   — sites in new but not original: shouldn't exist; flag if any

library(data.table)

args         <- commandArgs(trailingOnly = TRUE)
original_dir <- args[1]
new_dir      <- args[2]

chrs <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

# B-population founder sample names present in original RefAlt
# (used to exclude founder columns when comparing experimental pools)
founder_names <- c("AB8", "B1", "B2", "B3", "B4", "B5.RG", "B5", "B6", "B7")

results <- lapply(chrs, function(chr) {

  orig_file <- file.path(original_dir, paste0("RefAlt.", chr, ".txt"))
  new_file  <- file.path(new_dir,      paste0("RefAlt.", chr, ".txt"))

  if (!file.exists(orig_file)) { message("MISSING original: ", orig_file); return(NULL) }
  if (!file.exists(new_file))  { message("MISSING new: ",      new_file);  return(NULL) }

  message("\n=== ", chr, " ===")

  orig <- fread(orig_file, sep = "\t")
  new  <- fread(new_file,  sep = "\t")

  # Standardize column name separator (original uses spaces in data but tabs in header)
  setkey(orig, CHROM, POS)
  setkey(new,  CHROM, POS)

  orig_sites <- orig[, .(CHROM, POS)]
  new_sites  <- new[,  .(CHROM, POS)]

  n_orig   <- nrow(orig_sites)
  n_new    <- nrow(new_sites)
  shared   <- merge(orig_sites, new_sites, by = c("CHROM", "POS"))
  n_shared <- nrow(shared)

  only_orig <- fsetdiff(orig_sites, new_sites)
  only_new  <- fsetdiff(new_sites,  orig_sites)

  message(sprintf("  Original sites : %d", n_orig))
  message(sprintf("  New sites      : %d", n_new))
  message(sprintf("  Shared         : %d  (%.1f%% of original)", n_shared, 100*n_shared/n_orig))
  message(sprintf("  Only in original (dropped): %d", nrow(only_orig)))
  message(sprintf("  Only in new (gained):       %d", nrow(only_new)))

  if (nrow(only_new) > 0) {
    message("  *** WARNING: sites present in new but not original — unexpected ***")
    print(only_new[1:min(10, nrow(only_new))])
  }

  # --- Exact match check at shared sites ---
  # Find experimental pool columns (exclude CHROM, POS, and founder columns)
  is_founder_col <- function(col) {
    any(sapply(founder_names, function(f) grepl(paste0("_", f, "$"), col)))
  }
  pool_cols_orig <- names(orig)[!names(orig) %in% c("CHROM","POS") & !sapply(names(orig), is_founder_col)]
  pool_cols_new  <- names(new)[!names(new)   %in% c("CHROM","POS")]

  # Extract sample names from column headers (REF_<name> -> <name>)
  get_sample <- function(cols) unique(sub("^(REF|ALT)_", "", cols))
  orig_samples <- get_sample(pool_cols_orig)
  new_samples  <- get_sample(pool_cols_new)

  common_samples <- intersect(orig_samples, new_samples)
  message(sprintf("  Pool samples in original: %d", length(orig_samples)))
  message(sprintf("  Pool samples in new:      %d", length(new_samples)))
  message(sprintf("  Common pool samples:      %d", length(common_samples)))

  if (length(common_samples) == 0) {
    message("  *** No common sample names — check sample name convention ***")
    message("  Original examples: ", paste(head(orig_samples, 3), collapse=", "))
    message("  New examples:      ", paste(head(new_samples, 3), collapse=", "))
    return(NULL)
  }

  # Subset to shared sites
  orig_shared <- orig[shared, on = c("CHROM","POS")]
  new_shared  <- new[shared,  on = c("CHROM","POS")]

  # Check exact match per sample
  mismatch_summary <- rbindlist(lapply(common_samples, function(s) {
    ref_col <- paste0("REF_", s)
    alt_col <- paste0("ALT_", s)

    if (!ref_col %in% names(orig_shared) || !ref_col %in% names(new_shared)) return(NULL)

    ref_match <- all(orig_shared[[ref_col]] == new_shared[[ref_col]])
    alt_match <- all(orig_shared[[alt_col]] == new_shared[[alt_col]])
    n_ref_diff <- sum(orig_shared[[ref_col]] != new_shared[[ref_col]])
    n_alt_diff <- sum(orig_shared[[alt_col]] != new_shared[[alt_col]])

    data.table(sample=s, ref_exact=ref_match, alt_exact=alt_match,
               n_ref_diff=n_ref_diff, n_alt_diff=n_alt_diff)
  }))

  n_mismatch <- sum(!mismatch_summary$ref_exact | !mismatch_summary$alt_exact)
  if (n_mismatch == 0) {
    message(sprintf("  EXACT MATCH: all %d common samples identical at %d shared sites",
                    length(common_samples), n_shared))
  } else {
    message(sprintf("  *** MISMATCH in %d / %d samples — investigate ***", n_mismatch, length(common_samples)))
    print(mismatch_summary[!ref_exact | !alt_exact])
  }

  # --- Characterize dropped sites ---
  if (nrow(only_orig) > 0) {
    dropped <- orig[only_orig, on = c("CHROM","POS")]
    # Compute total read depth across all pool REF columns at dropped sites
    ref_cols_orig <- paste0("REF_", common_samples)
    ref_cols_orig <- ref_cols_orig[ref_cols_orig %in% names(dropped)]
    alt_cols_orig <- paste0("ALT_", common_samples)
    alt_cols_orig <- alt_cols_orig[alt_cols_orig %in% names(dropped)]

    if (length(ref_cols_orig) > 0) {
      dropped[, total_ref := rowSums(.SD), .SDcols = ref_cols_orig]
      dropped[, total_alt := rowSums(.SD), .SDcols = alt_cols_orig]
      dropped[, pool_alt_freq := total_alt / (total_ref + total_alt + 1e-6)]

      message("\n  Dropped sites — pool alt frequency distribution:")
      print(summary(dropped$pool_alt_freq))
      message("  (Sites with very low pool alt freq suggest they were driven by pools hitting QUAL>59)")

      # Check founder columns at dropped sites
      founder_ref_cols <- names(dropped)[grepl("^REF_", names(dropped)) & sapply(names(dropped), is_founder_col)]
      founder_alt_cols <- names(dropped)[grepl("^ALT_", names(dropped)) & sapply(names(dropped), is_founder_col)]
      if (length(founder_ref_cols) > 0) {
        dropped[, founder_total_ref := rowSums(.SD), .SDcols = founder_ref_cols]
        dropped[, founder_total_alt := rowSums(.SD), .SDcols = founder_alt_cols]
        dropped[, founder_alt_freq := founder_total_alt / (founder_total_ref + founder_total_alt + 1e-6)]
        message("  Dropped sites — founder alt frequency distribution:")
        print(summary(dropped$founder_alt_freq))
        message("  (Sites with low founder alt freq are marginal in founders too)")
      }
    }
  }

  list(chr=chr, n_orig=n_orig, n_new=n_new, n_shared=n_shared,
       n_dropped=nrow(only_orig), n_gained=nrow(only_new),
       mismatch_summary=mismatch_summary)
})

# --- Summary table ---
message("\n\n=== SUMMARY ===")
summary_dt <- rbindlist(lapply(results, function(r) {
  if (is.null(r)) return(NULL)
  data.table(chr=r$chr, n_orig=r$n_orig, n_new=r$n_new, n_shared=r$n_shared,
             n_dropped=r$n_dropped, n_gained=r$n_gained,
             pct_retained=round(100*r$n_shared/r$n_orig, 1))
}))
print(summary_dt)
