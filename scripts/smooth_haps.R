#!/usr/bin/env Rscript
###############################################################################
# smooth_haps.R  —  Step 1 of the freqsmooth pipeline
#
# Reads R.haps.<chr>.out.rds. Masks unresolvable founder frequencies to NA
# (founders that cluster together at a window have arbitrary individual
# estimates — only their sum is constrained). Then fills gaps and smooths:
#
#   1. Mask: set unresolvable founder frequencies to NA (per-founder, per-window)
#   2. Fill gaps: for each NA gap in a founder's series, compute the mean of
#      ~smooth_half resolved positions on each flank, then linearly interpolate
#      across the gap between those mean anchors.  Averaging many positions
#      avoids anchoring on the barely-resolved boundary values (which just
#      barely passed the hclust distance cutoff).
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
#       --rfile helpfiles/proj/design.txt --smooth-kb 125 --sex mixed
#
# --sex is the sex of the pools in THIS scan (M, F or mixed; default mixed).
# It only affects chrX -- see the chrX dosage section below.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args   <- commandArgs(trailingOnly = TRUE)
parsed <- list(smooth_kb = 125L, sex = "mixed")
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--chr"       = { parsed$chr       <- args[i+1]; i <- i+2L },
    "--dir"       = { parsed$dir       <- args[i+1]; i <- i+2L },
    "--outdir"    = { parsed$outdir    <- args[i+1]; i <- i+2L },
    "--rfile"     = { parsed$rfile     <- args[i+1]; i <- i+2L },
    "--smooth-kb" = { parsed$smooth_kb <- as.integer(args[i+1]); i <- i+2L },
    "--sex"       = { parsed$sex       <- args[i+1]; i <- i+2L },
    stop(paste("Unknown argument:", args[i]))
  )
}

mychr     <- parsed$chr
smooth_kb <- parsed$smooth_kb
design.df <- read.table(parsed$rfile, header = TRUE)

# chrX dosage.  Num counts flies, and everything downstream works in
# chromosomes: 2*Num of them.  That is right on an autosome, where every fly
# carries two.  On chrX the pool carries fewer, and how many fewer depends on
# the sex of the flies -- a female has 2 X, a male 1, so an equal mix averages
# 1.5.  Scaling Num by half the X-per-fly count makes 2*Num come out right:
#
#   F      2 X per fly  ->  1.00  ->  2.0*Num chromosomes
#   M      1 X per fly  ->  0.50  ->  1.0*Num
#   mixed  1.5 per fly  ->  0.75  ->  1.5*Num
#
# A scan is one sex throughout: contrasting a male pool against a female one
# would confound the treatment effect with sex, so --sex is per scan rather
# than per sample.  Default is mixed, which is the fixed value this pipeline
# applied unconditionally before, so omitting --sex reproduces old results
# exactly -- but it is wrong for a single-sex experiment (XQTL2 #38).
SEX_XFACTOR <- c(F = 1.0, M = 0.5, mixed = 0.75)
sex <- parsed$sex
if (!sex %in% names(SEX_XFACTOR))
  stop("--sex must be one of M, F or mixed (got '", sex, "').\n",
       "  M     = all-male pools     (1 X per fly)\n",
       "  F     = all-female pools   (2 X per fly)\n",
       "  mixed = equal numbers of males and females (the default)")
xfactor <- if (mychr == "chrX") unname(SEX_XFACTOR[sex]) else 1.0

# running_mean() and fill_gaps() live in scan_functions.R so this script and
# the covariance repair reuse them (XQTL2 #40).
script_dir <- dirname(normalizePath(sub("--file=", "",
                grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "scan_functions.R"))

# Stage layout: scans go under Scans/<scan> (see reorganize_project.sh).
dirout        <- file.path(parsed$dir, "Scans", parsed$outdir)
fileout_rds   <- file.path(dirout, paste0(parsed$outdir, ".smooth.", mychr, ".rds"))
fileout_means <- file.path(dirout, paste0(parsed$outdir, ".meansBySample.", mychr, ".txt"))
dir.create(dirout, showWarnings = FALSE, recursive = TRUE)

# ── Load ──────────────────────────────────────────────────────────────────────
filein <- file.path(parsed$dir, "Haps", paste0("R.haps.", mychr, ".out.rds"))

# Same guard REFALT2haps.R carries. This is the more likely place to hit it:
# REFALT2haps is skipped whenever haplotypes already exist -- the normal case for
# an old project, and the one the README recommends -- so the guarded script is
# the one that gets bypassed and this is where the operator lands (XQTL2 #31).
if (!file.exists(filein) &&
    file.exists(file.path(parsed$dir, paste0("R.haps.", mychr, ".out.rds"))))
  stop("This project uses the old flat layout (R.haps in ", parsed$dir, "/). Run:\n",
       "  bash pipeline/scripts/reorganize_project.sh ", parsed$dir, "\nthen rerun.")

cat("Reading", filein, "\n")
xx1 <- readRDS(filein)

step_bp     <- as.integer(median(diff(xx1$pos), na.rm = TRUE))
smooth_half <- round(smooth_kb * 1000L / step_bp)

cat(sprintf("  %d windows | step %d bp | smooth_half %d windows (+/-%d kb)\n",
            nrow(xx1), step_bp, smooth_half, smooth_kb))

if (mychr == "chrX")
  cat(sprintf("  chrX: --sex %s, so Num is scaled by %g (%g X chromosomes per fly)\n",
              sex, xfactor, 2*xfactor))

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
  mutate(Num = xfactor * Num)

# Mask unresolvable founders: if >1 founder shares a group at a window,
# those founders' individual frequencies are arbitrary (only their sum is
# constrained by the least-squares fit).  Set to NA so the interpolation +
# smoothing pipeline below recovers sensible values from flanking windows
# where the founders ARE resolved.
#
# The group sum IS known from lsei — save it so the constrained rescaling
# step below can restore it after gap-filling.
freq_raw <- freq_raw %>%
  group_by(CHROM, pos, pool, group) %>%
  mutate(group_size = n(),
         group_sum  = sum(freq, na.rm = TRUE)) %>%
  ungroup()

# Reference table: combined frequency per group at TRT/REP level
group_info <- freq_raw %>%
  filter(group_size > 1L) %>%
  group_by(CHROM, pos, TRT, REP, founder) %>%
  summarize(hclust_group = first(group),
            group_sum    = mean(group_sum, na.rm = TRUE),
            .groups      = "drop")

freq_raw <- freq_raw %>%
  mutate(freq = if_else(group_size > 1L, NA_real_, freq))

n_masked <- sum(is.na(freq_raw$freq))
n_total  <- nrow(freq_raw)
cat(sprintf("  Masked %d / %d founder-window estimates (%.1f%%) as unresolvable\n",
            n_masked, n_total, 100 * n_masked / n_total))

# Three-step frequency recovery (order matters):
#
#   1. fill_gaps() — interpolate NA gaps from mean-anchored flanks.
#      pmax uses na.rm=FALSE so masked NAs survive to fill_gaps.
#
#   2. Constrained rescaling — rescale grouped founders in-place so
#      their sum matches the original lsei combined (group_sum).
#      Done as a single group_by mutate; solo founders are untouched.
#
#   3. running_mean() — smooth the rescaled series.
freq_smoothed <- freq_raw %>%
  group_by(CHROM, pos, TRT, REP, founder) %>%
  summarize(freq = mean(freq, na.rm = TRUE),
            Num  = mean(Num,  na.rm = TRUE), .groups = "drop") %>%
  arrange(CHROM, pos) %>%
  mutate(freq = if_else(is.nan(freq), NA_real_, freq),
         freq = pmax(freq, 0.0003, na.rm = FALSE)) %>%
  group_by(TRT, REP, founder) %>%
  mutate(freq = fill_gaps(freq, smooth_half)) %>%
  ungroup() %>%
  left_join(group_info, by = c("CHROM", "pos", "TRT", "REP", "founder")) %>%
  group_by(CHROM, pos, TRT, REP, hclust_group) %>%
  mutate(freq = if_else(
    !is.na(hclust_group),
    freq / sum(freq, na.rm = TRUE) * first(group_sum),
    freq)) %>%
  ungroup() %>%
  select(-hclust_group, -group_sum) %>%
  arrange(CHROM, pos) %>%
  group_by(TRT, REP, founder) %>%
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

# Mask unresolvable founders in the covariance, exactly as the frequencies were
# masked above (XQTL2 #40).  Where >1 founder shares a group, lsei identifies
# only the group sum.  It does not identify their individual variances either:
# it returns a perfectly anti-correlated block.  Observed at chrX:9205000 in
# AGE_SY20_M_no89 --
#
#   Var(B6) = 0.5514   Var(B7) = 0.5516   Cov = -0.5515   corr = -1.0000
#   Var(B6) + Var(B7) = 1.1030      Var(B6 + B7) = 0.000023
#
# 0.5516 is not a possible variance for a quantity bounded in [0,1], where the
# maximum is 0.25.  These entries are wrong, not merely large: the linearised
# lsei covariance ignores the constraints the estimator itself is subject to.
#
# The frequencies for these founders are masked, gap-filled from flanking
# windows and rescaled to the lsei group sum, so what every downstream step
# actually squares or contrasts is an IMPUTED frequency.  The uncertainty that
# belongs beside it is the uncertainty of that imputation, which is what the
# flanking windows carry -- not the unbounded uncertainty of a raw estimate that
# was discarded.  So the covariance gets the same three steps: mask, fill_gaps,
# smooth.  An entry is masked if either founder it refers to is unresolved,
# since a covariance with an unidentified founder is itself unidentified.
#
# Index order: fi/fj index founders in the order they appear in Names within a
# (window, pool), which is the order of the Err matrix rows.
unresolved_tbl <- xx1 %>%
  select(CHROM, pos, sample, Names, Groups) %>%
  unnest(c(sample, Names, Groups)) %>%
  rename(pool = sample) %>%
  unnest(c(Names, Groups)) %>%
  group_by(CHROM, pos, pool) %>%
  mutate(idx = row_number()) %>%
  group_by(CHROM, pos, pool, Groups) %>%
  mutate(unresolved = n() > 1L) %>%
  ungroup() %>%
  select(CHROM, pos, pool, idx, unresolved)

err_long <- err_unnested %>%
  select(-Err) %>%
  tidyr::uncount(nF2) %>%
  mutate(fi = rep(fi_tmpl, nrow(err_unnested)),
         fj = rep(fj_tmpl, nrow(err_unnested)),
         v  = as.vector(t(v_mat))) %>%
  left_join(unresolved_tbl %>% rename(fi = idx, un_i = unresolved),
            by = c("CHROM", "pos", "pool", "fi")) %>%
  left_join(unresolved_tbl %>% rename(fj = idx, un_j = unresolved),
            by = c("CHROM", "pos", "pool", "fj")) %>%
  mutate(v = if_else(coalesce(un_i, FALSE) | coalesce(un_j, FALSE), NA_real_, v))

n_masked_cov <- sum(is.na(err_long$v))
cat(sprintf("  Masked %d / %d covariance entries (%.1f%%) as unresolvable\n",
            n_masked_cov, nrow(err_long), 100 * n_masked_cov / nrow(err_long)))

err_smoothed <- err_long %>%
  group_by(CHROM, pos, TRT, REP, fi, fj) %>%
  # na.rm=FALSE: a masked entry must survive the pool average to reach
  # fill_gaps, the same reason the frequency path masks before pmax(na.rm=FALSE).
  summarize(v = mean(v), .groups = "drop") %>%
  arrange(CHROM, pos) %>%
  group_by(TRT, REP, fi, fj) %>%
  mutate(v = fill_gaps(v, smooth_half),
         v = running_mean(v, smooth_half)) %>%
  ungroup()

# ── Save smoothed data for steps 2 and 3 ─────────────────────────────────────
founder_names <- sort(unique(freq_smoothed$founder))

cat("Writing smoothed RDS:", fileout_rds, "\n")
saveRDS(
  list(freq          = freq_smoothed,
       err           = err_smoothed,
       founder_names = founder_names,
       # What --sex this scan assumed, so the smoothed file records it rather
       # than it living only in whatever command line happened to produce it.
       sex           = sex,
       xfactor       = xfactor),
  fileout_rds
)

# ── meansBySample from smoothed frequencies ───────────────────────────────────
cat("Writing meansBySample:", fileout_means, "\n")
freq_smoothed %>%
  select(chr = CHROM, pos, TRT, REP, founder, freq) %>%
  filter(!is.na(freq)) %>%
  write.table(fileout_means)

cat("Done.\n")
