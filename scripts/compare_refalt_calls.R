# compare_refalt_calls.R — compare two RefAlt callsets that are NOT expected to
# be identical (e.g. the current QUAL-filtered pipeline vs the proposed founder-
# catalog pipeline). Unlike compare_refalt.sh (which checks byte-identity of the
# tiled caller), this quantifies how two different callsets relate.
#
# For each chromosome it answers three questions in plain terms:
#   1. How many SNPs does each keep, and how much do the site lists overlap?
#   2. At the SNPs both keep, do the per-sample REF/ALT counts agree?
#   3. What do the sites only one side keeps look like (count/frequency-wise)?
#
# Usage:
#   Rscript compare_refalt_calls.R <dirA_current> <dirB_proposed> [chr1,chr2,...]
#
# Writes a summary to stdout and, per chromosome, two TSVs in dirB:
#   compare.<chr>.a_only.txt   sites only in A (current)
#   compare.<chr>.b_only.txt   sites only in B (proposed)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
dirA <- args[1]   # current (reference)
dirB <- args[2]   # proposed
chrs <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

# long format: one row per (POS, sample) with REF/ALT counts and frequency
to_long <- function(path) {
  df <- fread(path, header = TRUE)   # RefAlt.<chr>.txt: CHROM POS REF_<s> ALT_<s> ...
  m <- melt(df, id.vars = c("CHROM", "POS"), variable.name = "lab", value.name = "count")
  m[, refalt := substr(lab, 1, 3)]
  m[, sample := substr(lab, 5, nchar(as.character(lab)))]
  w <- dcast(m, CHROM + POS + sample ~ refalt, value.var = "count")
  setnames(w, c("REF", "ALT"), c("ref", "alt"), skip_absent = TRUE)
  w[, N := ref + alt]
  w[, freq := fifelse(N > 0, ref / N, NA_real_)]
  w[]
}

cat(sprintf("Comparing A(current)=%s  vs  B(proposed)=%s\n\n", dirA, dirB))

for (chr in chrs) {
  fa <- file.path(dirA, paste0("RefAlt.", chr, ".txt"))
  fb <- file.path(dirB, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(fa) || !file.exists(fb)) {
    cat(sprintf("%-6s  SKIP (missing %s)\n", chr,
                paste(c(if (!file.exists(fa)) "A", if (!file.exists(fb)) "B"), collapse = "+")))
    next
  }

  A <- to_long(fa); B <- to_long(fb)
  posA <- unique(A[, .(POS)]); posB <- unique(B[, .(POS)])
  nA <- nrow(posA); nB <- nrow(posB)
  shared <- fintersect(posA, posB)
  nShared <- nrow(shared)

  # (2) agreement at shared sites: join per (POS, sample), compare counts
  j <- merge(A, B, by = c("CHROM", "POS", "sample"), suffixes = c(".a", ".b"))
  identical_counts <- j[, mean(ref.a == ref.b & alt.a == alt.b)]
  # frequency difference where both have coverage
  jf <- j[!is.na(freq.a) & !is.na(freq.b)]
  mad_freq <- if (nrow(jf)) jf[, mean(abs(freq.a - freq.b))] else NA_real_

  cat(sprintf("%-6s  A=%d  B=%d  shared=%d (%.1f%% of A, %.1f%% of B)\n",
              chr, nA, nB, nShared, 100 * nShared / nA, 100 * nShared / nB))
  cat(sprintf("        at shared (POS,sample): %.2f%% identical counts; mean |freq diff| = %s\n",
              100 * identical_counts,
              ifelse(is.na(mad_freq), "NA", sprintf("%.4f", mad_freq))))

  # (3) dump the sites unique to each side for inspection
  aOnly <- fsetdiff(posA, posB); bOnly <- fsetdiff(posB, posA)
  fwrite(A[aOnly, on = "POS"], file.path(dirB, sprintf("compare.%s.a_only.txt", chr)), sep = "\t")
  fwrite(B[bOnly, on = "POS"], file.path(dirB, sprintf("compare.%s.b_only.txt", chr)), sep = "\t")
  cat(sprintf("        A-only=%d  B-only=%d  (rows written to %s/compare.%s.*_only.txt)\n\n",
              nrow(aOnly), nrow(bOnly), dirB, chr))
}
