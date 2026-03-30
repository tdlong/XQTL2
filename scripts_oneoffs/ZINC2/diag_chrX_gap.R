#!/usr/bin/env Rscript
###############################################################################
# diag_chrX_gap.R  —  Diagnose the chrX 20.5-22.6 Mb coverage gap
#
# Reads R.haps.chrX.out.rds and for each window in the 19-23 Mb region shows:
#   - Which founders are unresolvable (group_size > 1)
#   - Whether the gap is wider than the ±250 kb smoothing kernel
#     (if so, running_mean returns NA and the window is dropped from output)
#
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

HAPS_FILE <- "process/ZINC2/R.haps.chrX.out.rds"
DESIGN    <- "helpfiles/ZINC2/Zinc2.test.F.txt"
SMOOTH_KB <- 250
STEP_BP   <- 5000
SMOOTH_HALF_WINDOWS <- round(SMOOTH_KB * 1000 / STEP_BP)   # 50 windows

cat("Reading", HAPS_FILE, "\n")
xx1 <- readRDS(HAPS_FILE)
design.df <- read.table(DESIGN, header = TRUE)

cat(sprintf("  %d windows on chrX\n\n", nrow(xx1)))

# Unnest one representative pool (use first control sample) to get group structure
# group_size > 1 means founders are indistinguishable (arbitrary individual freqs)
one_pool <- design.df$bam[design.df$TRT == "C"][1]
cat("Using pool:", one_pool, "to characterise group structure\n\n")

groups_df <- xx1 %>%
  select(CHROM, pos, sample, Names, Groups) %>%
  unnest(c(sample, Names, Groups)) %>%
  unnest(c(Names, Groups)) %>%
  rename(pool = sample, founder = Names, group = Groups) %>%
  filter(pool == one_pool) %>%
  group_by(CHROM, pos, group) %>%
  mutate(group_size = n()) %>%
  ungroup()

# Per-window: number of masked founders and their names
win_mask <- groups_df %>%
  group_by(pos) %>%
  summarise(
    n_founders    = n(),
    n_masked      = sum(group_size > 1),
    masked_founders = paste(sort(founder[group_size > 1]), collapse = ","),
    all_masked    = all(group_size > 1),
    .groups = "drop"
  ) %>%
  arrange(pos) %>%
  mutate(pos_mb = pos / 1e6)

# Show 18-23 Mb region
region <- win_mask %>% filter(pos_mb >= 18, pos_mb <= 23)
cat("=== Masking per window, 18-23 Mb ===\n")
cat(sprintf("%-10s  %-5s  %-8s  %-11s  %s\n",
            "pos_mb", "n_fnd", "n_masked", "all_masked", "masked_founders"))
for (i in seq_len(nrow(region))) {
  r <- region[i, ]
  cat(sprintf("%-10.3f  %-5d  %-8d  %-11s  %s\n",
              r$pos_mb, r$n_founders, r$n_masked,
              as.character(r$all_masked), r$masked_founders))
}

# Find contiguous all-masked stretch
all_masked_windows <- win_mask %>% filter(all_masked)
cat(sprintf("\n=== Contiguous all-masked windows ===\n"))
if (nrow(all_masked_windows) == 0) {
  cat("None found.\n")
} else {
  cat(sprintf("  All-masked windows: %d\n", nrow(all_masked_windows)))
  cat(sprintf("  First: %.3f Mb\n", min(all_masked_windows$pos_mb)))
  cat(sprintf("  Last:  %.3f Mb\n", max(all_masked_windows$pos_mb)))
  cat(sprintf("  Width: %.2f Mb\n",
              (max(all_masked_windows$pos) - min(all_masked_windows$pos)) / 1e6))
  cat(sprintf("  Smoothing kernel: +/-%d kb = +/-%d windows\n",
              SMOOTH_KB, SMOOTH_HALF_WINDOWS))
  gap_half_kb <- (max(all_masked_windows$pos) - min(all_masked_windows$pos)) / 2000
  cat(sprintf("  Half-gap: %.1f kb vs kernel: %d kb\n", gap_half_kb, SMOOTH_KB))
  cat(sprintf("  Gap exceeds kernel? %s\n",
              if (gap_half_kb > SMOOTH_KB) "YES — interpolation fails, windows dropped" else "no"))
}

cat("\nDone.\n")
