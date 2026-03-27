#!/usr/bin/env Rscript
###############################################################################
# smooth_haps.R  —  Step 1 of the freqsmooth pipeline
#
# Reads R.haps.<chr>.out.rds. Masks unresolvable founder frequencies to NA
# (founders that cluster together at a window have arbitrary individual
# estimates — only their sum is constrained). Then fills gaps and smooths:
#
#   1. Mask: set unresolvable founder frequencies to NA (per-founder, per-window)
#   2. Fill gaps: for each NA gap in a founder's series, fit a linear trend
#      from ~smooth_half resolved positions on each flank, extrapolate to the
#      gap edges, then linearly interpolate across the gap.  This avoids
#      anchoring on the barely-resolved positions right at the gap boundary
#      (which just barely passed the hclust distance cutoff).
#   3. Smooth: apply a running mean (+/- smooth_half windows) to the filled
#      series.
#
# Saves:
#   <outdir>/<scan>.smooth.<chr>.rds       -- smoothed data for steps 2 & 3
#   <outdir>/<scan>.meansBySample.<chr>.txt -- smoothed per-founder frequencies
#
# Usage:
#   Rscript scripts/smooth_haps.R \
#       --chr chrX --dir process/proj --outdir SCAN_NAME \
#       --rfile helpfiles/proj/design.txt --smooth-kb 125
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args   <- commandArgs(trailingOnly = TRUE)
parsed <- list(smooth_kb = 125L)
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--chr"       = { parsed$chr       <- args[i+1]; i <- i+2L },
    "--dir"       = { parsed$dir       <- args[i+1]; i <- i+2L },
    "--outdir"    = { parsed$outdir    <- args[i+1]; i <- i+2L },
    "--rfile"     = { parsed$rfile     <- args[i+1]; i <- i+2L },
    "--smooth-kb" = { parsed$smooth_kb <- as.integer(args[i+1]); i <- i+2L },
    stop(paste("Unknown argument:", args[i]))
  )
}

mychr     <- parsed$chr
smooth_kb <- parsed$smooth_kb
design.df <- read.table(parsed$rfile, header = TRUE)

dirout        <- file.path(parsed$dir, parsed$outdir)
fileout_rds   <- file.path(dirout, paste0(parsed$outdir, ".smooth.", mychr, ".rds"))
fileout_means <- file.path(dirout, paste0(parsed$outdir, ".meansBySample.", mychr, ".txt"))
dir.create(dirout, showWarnings = FALSE, recursive = TRUE)

# ── Running mean (edge-aware, O(n), NA-safe) ──────────────────────────────────
running_mean <- function(x, h) {
  n  <- length(x); if (n == 0L) return(numeric(0))
  ok <- !is.na(x)
  xc <- replace(x, !ok, 0)
  cs <- c(0, cumsum(xc)); cn <- c(0, cumsum(as.numeric(ok)))
  lo <- pmax(1L, seq_len(n) - h); hi <- pmin(n, seq_len(n) + h)
  tot <- cs[hi+1L] - cs[lo]; cnt <- cn[hi+1L] - cn[lo]
  ifelse(cnt > 0L, tot/cnt, NA_real_)
}

# ── Gap filler (trend-based interpolation) ───────────────────────────────────
# For each contiguous run of NAs ("gap"), fit a linear trend from up to h
# valid positions on each flank, extrapolate to the gap edges, then linearly
# interpolate across the gap.
#
# Why not just use the values right at the gap boundary?  The haplotype
# estimator resolves founders by cutting a distance tree (hclust + cutree).
# Positions right at the gap edge just barely passed the cutoff — their
# frequency estimates are only marginally better than the unresolved ones
# inside the gap.  Fitting a trend from h flanking positions gives robust
# anchor values driven by the well-resolved interior, not the noisy boundary.
#
# Leading/trailing NAs (no flank on one side) get a flat extrapolation from
# the available flank's trend.  If no valid data exists at all, NAs remain.

fill_gaps <- function(x, h) {
  n  <- length(x)
  if (n == 0L || !anyNA(x)) return(x)

  ok  <- !is.na(x)
  if (!any(ok)) return(x)                   # all NA — nothing to anchor on

  # Identify contiguous NA runs (gaps)
  rle_na  <- rle(!ok)
  ends    <- cumsum(rle_na$lengths)
  starts  <- ends - rle_na$lengths + 1L

  for (g in which(rle_na$values)) {
    ga <- starts[g]                          # first NA in this gap
    gb <- ends[g]                            # last  NA in this gap

    # Left flank: up to h valid positions before the gap
    left_idx <- which(ok & seq_along(x) < ga)
    if (length(left_idx) > 0L) {
      left_idx <- tail(left_idx, h)
      if (length(left_idx) >= 2L) {
        fit_L <- lm(x[left_idx] ~ left_idx)
        anchor_L <- predict(fit_L, newdata = data.frame(left_idx = ga))
      } else {
        anchor_L <- x[left_idx]              # single point — use as-is
      }
    } else {
      anchor_L <- NULL
    }

    # Right flank: up to h valid positions after the gap
    right_idx <- which(ok & seq_along(x) > gb)
    if (length(right_idx) > 0L) {
      right_idx <- head(right_idx, h)
      if (length(right_idx) >= 2L) {
        fit_R <- lm(x[right_idx] ~ right_idx)
        anchor_R <- predict(fit_R, newdata = data.frame(right_idx = gb))
      } else {
        anchor_R <- x[right_idx]
      }
    } else {
      anchor_R <- NULL
    }

    # Fill the gap
    gap_len <- gb - ga + 1L
    if (!is.null(anchor_L) && !is.null(anchor_R)) {
      # Both flanks: linear interpolation between the two trend-based anchors
      x[ga:gb] <- seq(anchor_L, anchor_R, length.out = gap_len)
    } else if (!is.null(anchor_L)) {
      # Leading edge only: flat fill from left trend
      x[ga:gb] <- anchor_L
    } else if (!is.null(anchor_R)) {
      # Trailing edge only: flat fill from right trend
      x[ga:gb] <- anchor_R
    }
    # else: no flanks at all, leave as NA
  }
  x
}

# ── Load ──────────────────────────────────────────────────────────────────────
filein <- file.path(parsed$dir, paste0("R.haps.", mychr, ".out.rds"))
cat("Reading", filein, "\n")
xx1 <- readRDS(filein)

step_bp     <- as.integer(median(diff(xx1$pos), na.rm = TRUE))
smooth_half <- round(smooth_kb * 1000L / step_bp)
sexlink     <- if (mychr == "chrX") 0.75 else 1.0

cat(sprintf("  %d windows | step %d bp | smooth_half %d windows (+/-%d kb)\n",
            nrow(xx1), step_bp, smooth_half, smooth_kb))

options(dplyr.summarise.inform = FALSE)

# ── Smooth haplotype frequencies ──────────────────────────────────────────────
# Unnest: one row per (window, pool, founder)
# Average within (window, TRT, REP, founder) to collapse any technical reps
# Then group_by(TRT, REP, founder) and mutate running_mean across windows
cat("Smoothing frequencies...\n")

freq_raw <- xx1 %>%
  select(CHROM, pos, sample, Haps, Names, Groups) %>%
  unnest(c(sample, Haps, Names, Groups)) %>%
  unnest(c(Haps, Names, Groups)) %>%
  rename(pool = sample, freq = Haps, founder = Names, group = Groups) %>%
  left_join(design.df, by = c("pool" = "bam")) %>%
  filter(!is.na(TRT)) %>%
  mutate(Num = sexlink * Num)

# Mask unresolvable founders: if >1 founder shares a group at a window,
# those founders' individual frequencies are arbitrary (only their sum is
# constrained by the least-squares fit).  Set to NA so the interpolation +
# smoothing pipeline below recovers sensible values from flanking windows
# where the founders ARE resolved.
freq_raw <- freq_raw %>%
  group_by(CHROM, pos, pool, group) %>%
  mutate(group_size = n()) %>%
  ungroup() %>%
  mutate(freq = if_else(group_size > 1L, NA_real_, freq))

n_masked <- sum(is.na(freq_raw$freq))
n_total  <- nrow(freq_raw)
cat(sprintf("  Masked %d / %d founder-window estimates (%.1f%%) as unresolvable\n",
            n_masked, n_total, 100 * n_masked / n_total))

# Two-step frequency recovery (order matters — fill gaps THEN smooth):
#
#   1. fill_gaps() — for each founder's NA gaps, fit a linear trend from ~h
#      resolved positions on each flank, extrapolate to the gap edges, then
#      linearly interpolate across the gap.  This uses the trend from the
#      well-resolved interior (not just the barely-resolved boundary positions)
#      to anchor the fill.
#
#   2. running_mean() — smooth the now-complete series with a +/- smooth_half
#      window.  Because gaps are already filled, the smoother sees a continuous
#      series and produces uniform-quality output everywhere.
freq_smoothed <- freq_raw %>%
  group_by(CHROM, pos, TRT, REP, founder) %>%
  summarize(freq = mean(freq, na.rm = TRUE),
            Num  = mean(Num,  na.rm = TRUE), .groups = "drop") %>%
  arrange(CHROM, pos) %>%
  group_by(TRT, REP, founder) %>%
  mutate(freq = fill_gaps(freq, smooth_half)) %>%
  mutate(freq = running_mean(freq, smooth_half)) %>%
  ungroup()

# ── Smooth reconstruction covariance matrices ─────────────────────────────────
# Unnest: one row per (window, pool), then expand each matrix to (fi, fj, v) rows
# Average within (window, TRT, REP, fi, fj) then smooth across windows
cat("Smoothing covariance matrices...\n")

err_unnested <- xx1 %>%
  select(CHROM, pos, sample, Err) %>%
  unnest(c(sample, Err)) %>%
  rename(pool = sample) %>%
  left_join(design.df %>% select(bam, TRT, REP), by = c("pool" = "bam")) %>%
  filter(!is.na(TRT))

# Vectorized covariance expansion — column-major order matches as.vector()
nF_cov  <- nrow(as.matrix(err_unnested$Err[[1]]))
nF2     <- nF_cov^2L
fi_tmpl <- rep(seq_len(nF_cov), nF_cov)               # row indices
fj_tmpl <- rep(seq_len(nF_cov), each = nF_cov)        # col indices
v_mat   <- do.call(rbind, lapply(err_unnested$Err,     # n_rows x nF2
             function(m) as.vector(as.matrix(m))))

err_smoothed <- err_unnested %>%
  select(-Err) %>%
  tidyr::uncount(nF2) %>%
  mutate(fi = rep(fi_tmpl, nrow(err_unnested)),
         fj = rep(fj_tmpl, nrow(err_unnested)),
         v  = as.vector(t(v_mat))) %>%
  group_by(CHROM, pos, TRT, REP, fi, fj) %>%
  summarize(v = mean(v, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHROM, pos) %>%
  group_by(TRT, REP, fi, fj) %>%
  mutate(v = running_mean(v, smooth_half)) %>%
  ungroup()

# ── Save smoothed data for steps 2 and 3 ─────────────────────────────────────
founder_names <- sort(unique(freq_smoothed$founder))
nrepl         <- n_distinct(freq_smoothed$REP)

cat("Writing smoothed RDS:", fileout_rds, "\n")
saveRDS(
  list(freq          = freq_smoothed,
       err           = err_smoothed,
       founder_names = founder_names,
       nrepl         = nrepl),
  fileout_rds
)

# ── meansBySample from smoothed frequencies ───────────────────────────────────
cat("Writing meansBySample:", fileout_means, "\n")
freq_smoothed %>%
  select(chr = CHROM, pos, TRT, REP, founder, freq) %>%
  filter(!is.na(freq)) %>%
  write.table(fileout_means)

cat("Done.\n")
