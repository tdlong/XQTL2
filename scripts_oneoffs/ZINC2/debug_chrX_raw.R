#!/usr/bin/env Rscript
###############################################################################
# debug_chrX_raw.R — Check what fill_gaps sees and produces
#
# Replays the smooth_haps.R pipeline for chrX, dumping intermediate values
# for B2 and B4 at the base (>20Mb):
#   1. Raw freq after masking (before fill_gaps)
#   2. After fill_gaps
#   3. After running_mean
#
# Based on chrX_freq_check.R (which works).
# Run from: /dfs7/adl/tdlong/fly_pool/XQTL2
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

RDSFILE <- "process/ZINC2/R.haps.chrX.out.rds"
DESIGN  <- "helpfiles/ZINC2/Zinc2.test.F.txt"
OUTDIR  <- "scripts_oneoffs/ZINC2/debug_chrX_results"
OUTTEXT <- file.path(OUTDIR, "raw_diag.txt")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
sink(OUTTEXT, split = TRUE)
options(width = 200)

smooth_kb <- 250L
design.df <- read.table(DESIGN, header = TRUE)

# ── Replicate smooth_haps.R functions ────────────────────────────────────────
running_mean <- function(x, h) {
  n  <- length(x); if (n == 0L) return(numeric(0))
  ok <- !is.na(x)
  xc <- replace(x, !ok, 0)
  cs <- c(0, cumsum(xc)); cn <- c(0, cumsum(as.numeric(ok)))
  lo <- pmax(1L, seq_len(n) - h); hi <- pmin(n, seq_len(n) + h)
  tot <- cs[hi+1L] - cs[lo]; cnt <- cn[hi+1L] - cn[lo]
  ifelse(cnt > 0L, tot/cnt, NA_real_)
}

fill_gaps <- function(x, h) {
  n  <- length(x)
  if (n == 0L || !anyNA(x)) return(x)
  ok  <- !is.na(x)
  if (!any(ok)) return(x)
  rle_na  <- rle(!ok)
  ends    <- cumsum(rle_na$lengths)
  starts  <- ends - rle_na$lengths + 1L
  for (g in which(rle_na$values)) {
    ga <- starts[g]; gb <- ends[g]
    left_idx <- which(ok & seq_along(x) < ga)
    if (length(left_idx) > 0L) {
      anchor_L <- mean(x[tail(left_idx, h)])
    } else { anchor_L <- NULL }
    right_idx <- which(ok & seq_along(x) > gb)
    if (length(right_idx) > 0L) {
      anchor_R <- mean(x[head(right_idx, h)])
    } else { anchor_R <- NULL }
    gap_len <- gb - ga + 1L
    if (!is.null(anchor_L) && !is.null(anchor_R)) {
      x[ga:gb] <- seq(anchor_L, anchor_R, length.out = gap_len)
    } else if (!is.null(anchor_L)) {
      x[ga:gb] <- anchor_L
    } else if (!is.null(anchor_R)) {
      x[ga:gb] <- anchor_R
    }
  }
  x
}

# ── Load and unnest ─────────────────────────────────────────────────────────
cat("Reading", RDSFILE, "\n")
xx1 <- readRDS(RDSFILE)
step_bp     <- as.integer(median(diff(xx1$pos), na.rm = TRUE))
smooth_half <- round(smooth_kb * 1000L / step_bp)
cat(sprintf("  %d windows | step %d bp | smooth_half %d\n", nrow(xx1), step_bp, smooth_half))

freq_raw <- xx1 %>%
  select(CHROM, pos, sample, Haps, Names, Groups) %>%
  unnest(c(sample, Haps, Names, Groups)) %>%
  unnest(c(Haps, Names, Groups)) %>%
  rename(pool = sample, freq = Haps, founder = Names, group = Groups) %>%
  left_join(design.df, by = c("pool" = "bam")) %>%
  filter(!is.na(TRT)) %>%
  mutate(Num = Num)

# ── Check: what are the raw (pre-mask) values for B2/B4 at the base? ────────
cat("\n=== RAW freq (before masking) for B2/B4, pos > 20Mb, REP=1, TRT=C ===\n")
raw_check <- freq_raw %>%
  filter(founder %in% c("B2", "B4"), pos > 20e6, REP == 1, TRT == "C") %>%
  group_by(founder, pos) %>%
  summarize(n_pools = n(), mean_freq = mean(freq), min_freq = min(freq),
            max_freq = max(freq), .groups = "drop") %>%
  filter(pos > 21e6, pos < 22.5e6) %>%
  arrange(founder, pos)
print(as.data.frame(raw_check), row.names = FALSE)

# ── Mask ─────────────────────────────────────────────────────────────────────
freq_raw <- freq_raw %>%
  group_by(CHROM, pos, pool, group) %>%
  mutate(group_size = n()) %>%
  ungroup() %>%
  mutate(freq = if_else(group_size > 1L, NA_real_, freq))

cat("\n=== After masking, B2/B4 at base, REP=1, TRT=C ===\n")
masked_check <- freq_raw %>%
  filter(founder %in% c("B2", "B4"), pos > 20e6, REP == 1, TRT == "C") %>%
  group_by(founder, pos) %>%
  summarize(n_valid = sum(!is.na(freq)), mean_freq = mean(freq, na.rm = TRUE),
            .groups = "drop") %>%
  filter(pos > 21e6, pos < 22.5e6) %>%
  arrange(founder, pos)
print(as.data.frame(masked_check), row.names = FALSE)

# ── Summarize (collapse pools) then fill gaps ────────────────────────────────
freq_summ <- freq_raw %>%
  group_by(CHROM, pos, TRT, REP, founder) %>%
  summarize(freq = mean(freq, na.rm = TRUE),
            Num  = mean(Num, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHROM, pos)

cat("\n=== After summarize (pool mean), B2/B4 at base, REP=1, TRT=C ===\n")
summ_check <- freq_summ %>%
  filter(founder %in% c("B2", "B4"), pos > 21e6, pos < 22.5e6,
         REP == 1, TRT == "C") %>%
  arrange(founder, pos)
print(as.data.frame(summ_check), row.names = FALSE)

# Apply fill_gaps
freq_filled <- freq_summ %>%
  group_by(TRT, REP, founder) %>%
  mutate(freq_filled = fill_gaps(freq, smooth_half)) %>%
  ungroup()

cat("\n=== After fill_gaps, B2/B4 at base, REP=1, TRT=C ===\n")
fill_check <- freq_filled %>%
  filter(founder %in% c("B2", "B4"), pos > 21e6, pos < 22.5e6,
         REP == 1, TRT == "C") %>%
  select(founder, pos, freq_before = freq, freq_after = freq_filled) %>%
  arrange(founder, pos)
print(as.data.frame(fill_check), row.names = FALSE)

cat("\nDone.\n")
sink()

# Auto-push
system("cd /dfs7/adl/tdlong/fly_pool/XQTL2 && git add scripts_oneoffs/ZINC2/debug_chrX_results/raw_diag.txt && git commit -m 'debug raw_diag results' && git push dev HEAD:main 2>&1")
