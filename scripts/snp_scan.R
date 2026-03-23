#!/usr/bin/env Rscript
###############################################################################
# snp_scan.R  —  Step 3 of the freqsmooth pipeline
#
# For each SNP, imputes ALT frequency in each pool as:
#
#   f_ALT(pool, SNP) = h(pool, window)  ·  s(SNP)
#
# where h is the smoothed 8-founder frequency vector and s is the vector of
# founder ALT states from FREQ_SNPs.cM.txt.gz.
#
# SNPs are grouped by their nearest haplotype window. Within each window group,
# the frequency computation is a single matrix multiply (vectorized over all
# SNPs in the group), and the propagated covariance is computed likewise.
# The stabilized Wald test (df=1) follows.
#
# Usage:
#   Rscript scripts/snp_scan.R \
#       --chr chrX --dir process/proj/SCAN_NAME --outdir SCAN_NAME \
#       --rfile helpfiles/proj/design.txt \
#       --snp-table FREQ_SNPs.cM.txt.gz \
#       --founders A1,A2,A3,A4,A5,A6,A7,AB8
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

RIDGE_FRACTION <- 0.01
AF_CUTOFF      <- 0.01

# ── Arguments ─────────────────────────────────────────────────────────────────
args   <- commandArgs(trailingOnly = TRUE)
parsed <- list()
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--chr"       = { parsed$chr       <- args[i+1];                              i <- i+2L },
    "--dir"       = { parsed$dir       <- args[i+1];                              i <- i+2L },
    "--outdir"    = { parsed$outdir    <- args[i+1];                              i <- i+2L },
    "--rfile"     = { parsed$rfile     <- args[i+1];                              i <- i+2L },
    "--snp-table" = { parsed$snp_table <- args[i+1];                              i <- i+2L },
    "--founders"  = { parsed$founders  <- strsplit(args[i+1], ",")[[1]];          i <- i+2L },
    stop(paste("Unknown argument:", args[i]))
  )
}

mychr    <- parsed$chr
founders <- parsed$founders
design.df <- read.table(parsed$rfile, header = TRUE)
source("scripts/scan_functions.R")

filein        <- file.path(parsed$dir, paste0(parsed$outdir, ".smooth.", mychr, ".rds"))
fileout       <- file.path(parsed$dir, paste0(parsed$outdir, ".snp_scan.", mychr, ".txt"))
fileout_means <- file.path(parsed$dir, paste0(parsed$outdir, ".snp_meansBySample.", mychr, ".txt"))

options(dplyr.summarise.inform = FALSE)

# ── Load smoothed haplotype data ──────────────────────────────────────────────
cat("Reading", filein, "\n")
sm            <- readRDS(filein)
freq_smoothed <- sm$freq
err_smoothed  <- sm$err
nrepl         <- sm$nrepl
nF            <- length(founders)

# Window positions for nearest-join
win_positions <- freq_smoothed %>%
  distinct(CHROM, pos) %>%
  filter(CHROM == mychr) %>%
  arrange(pos) %>%
  pull(pos)
W <- length(win_positions)

# ── Load SNP table for this chromosome ───────────────────────────────────────
cat("Reading SNP table for", mychr, "...\n")
snp_df <- read.table(parsed$snp_table, header = TRUE) %>%
  filter(CHROM == mychr)
cat(sprintf("  %d SNPs\n", nrow(snp_df)))

if (nrow(snp_df) == 0) {
  tibble(chr=character(), pos=integer(), Wald_log10p=numeric(),
         cM=numeric(), n_informative_founders=integer()) %>%
    write.table(fileout)
  tibble(chr=character(), pos=integer(), TRT=character(),
         REP=integer(), F_alt=numeric(), cM=numeric()) %>%
    write.table(fileout_means)
  quit(save="no")
}

# Founder state matrix: n_SNPs x nF
S_mat <- as.matrix(snp_df[, founders])

# Assign each SNP to its nearest haplotype window (vectorized via findInterval)
mid_pts  <- (win_positions[-1] + win_positions[-W]) / 2
win_idx  <- findInterval(snp_df$POS, mid_pts) + 1L
win_idx  <- pmin(pmax(win_idx, 1L), W)
snp_df$win_idx <- win_idx

# ── Pre-build per-window frequency and covariance objects ────────────────────
# For each window: H_C[r, f] and H_Z[r, f] (nrepl x nF freq matrices)
#                  Err_C[f,f,r] and Err_Z[f,f,r] (nF x nF x nrepl cov arrays)
cat("Pre-building per-window arrays...\n")

freq_chr <- freq_smoothed %>% filter(CHROM == mychr)
err_chr  <- err_smoothed  %>% filter(CHROM == mychr)

# Nest by position for fast lookup
freq_by_pos <- freq_chr %>% group_by(pos) %>% group_split() %>%
  setNames(win_positions)
err_by_pos  <- err_chr  %>% group_by(pos) %>% group_split() %>%
  setNames(win_positions)

build_freq_mat <- function(w_freq, trt) {
  w_freq %>% filter(TRT == trt) %>% arrange(REP) %>%
    pivot_wider(names_from = founder, values_from = freq) %>%
    select(all_of(founders)) %>% as.matrix()
}

build_err_arr <- function(w_err, trt) {
  arr <- array(NA_real_, c(nF, nF, nrepl))
  for (r in seq_len(nrepl)) {
    vals <- w_err %>% filter(TRT == trt, REP == r) %>% arrange(fi, fj) %>% pull(v)
    if (length(vals) == nF^2) arr[,,r] <- matrix(vals, nF, nF)
  }
  arr
}

# Cache N1/N2 (constant across windows for a given REP)
N1 <- freq_chr %>% filter(TRT=="C", CHROM==mychr) %>%
  group_by(REP) %>% slice(1) %>% arrange(REP) %>% pull(Num)
N2 <- freq_chr %>% filter(TRT=="Z", CHROM==mychr) %>%
  group_by(REP) %>% slice(1) %>% arrange(REP) %>% pull(Num)

# ── SNP scan: process one window group at a time ─────────────────────────────
# All SNPs in the same window share the same H and Err, so the frequency
# computation is a matrix multiply over all SNPs at once.
cat(sprintf("Running SNP scan (%d SNPs across %d windows)...\n",
            nrow(snp_df), W))

all_results <- snp_df %>%
  group_by(win_idx) %>%
  group_map(function(snps, key) {

    w_pos  <- win_positions[key$win_idx]
    w_freq <- freq_by_pos[[as.character(w_pos)]]
    w_err  <- err_by_pos [[as.character(w_pos)]]
    if (is.null(w_freq) || is.null(w_err)) return(NULL)

    H_C <- build_freq_mat(w_freq, "C")   # nrepl x nF
    H_Z <- build_freq_mat(w_freq, "Z")
    if (anyNA(H_C) || anyNA(H_Z)) return(NULL)

    E_C <- build_err_arr(w_err, "C")     # nF x nF x nrepl
    E_Z <- build_err_arr(w_err, "Z")
    if (anyNA(E_C) || anyNA(E_Z)) return(NULL)

    # Founder state matrix for this group of SNPs: n_snps x nF
    S <- as.matrix(snps[, founders])

    # Filter to SNPs polymorphic across founders
    n_alt  <- rowSums(S > 0.5)
    n_ref  <- rowSums(S < 0.5)
    inform <- n_alt > 0 & n_ref > 0
    if (!any(inform)) return(NULL)
    S      <- S[inform, , drop = FALSE]
    snps_i <- snps[inform, ]
    n_s    <- nrow(S)

    # ── Vectorized SNP freq: all SNPs × all replicates at once ───────────────
    # F_alt[s, r] = S[s,] · H[r,]  →  S %*% t(H)   (n_s × nrepl)
    F_alt_C <- S %*% t(H_C)
    F_alt_Z <- S %*% t(H_Z)

    # Broadcast N1/N2 to n_s × nrepl
    N1_mat <- matrix(N1, n_s, nrepl, byrow = TRUE)
    N2_mat <- matrix(N2, n_s, nrepl, byrow = TRUE)

    # Pool-mean ALT freq: n_s × nrepl
    pm_alt <- (N1_mat * F_alt_C + N2_mat * F_alt_Z) / (N1_mat + N2_mat)

    # Covariance propagation — loop is over nrepl only, not over SNPs
    # var_alt[s,r] = s · E[,,r] · s  =  rowSums( (S %*% E[,,r]) * S )
    Sm1       <- 1 - S
    var_alt_C <- matrix(0, n_s, nrepl)
    var_alt_Z <- matrix(0, n_s, nrepl)
    var_ref_C <- matrix(0, n_s, nrepl)
    var_ref_Z <- matrix(0, n_s, nrepl)
    for (r in seq_len(nrepl)) {
      var_alt_C[,r] <- rowSums((S   %*% E_C[,,r]) * S)
      var_alt_Z[,r] <- rowSums((S   %*% E_Z[,,r]) * S)
      var_ref_C[,r] <- rowSums((Sm1 %*% E_C[,,r]) * Sm1)
      var_ref_Z[,r] <- rowSums((Sm1 %*% E_Z[,,r]) * Sm1)
    }

    # sum(diag(mn.covmat(pm, 2*N))) for biallelic = pm*(1-pm)/N
    sdcm1 <- pm_alt * (1 - pm_alt) / N1_mat
    sdcm2 <- pm_alt * (1 - pm_alt) / N2_mat

    # Effective N: n_s × nrepl
    N1e <- 2 * sdcm1 * N1_mat / (sdcm1 + var_alt_C + var_ref_C)
    N2e <- 2 * sdcm2 * N2_mat / (sdcm2 + var_alt_Z + var_ref_Z)

    # Pooled frequencies: n_s vectors
    sum_N1e <- rowSums(N1e)
    sum_N2e <- rowSums(N2e)
    p1p <- rowSums(N1e * F_alt_C) / sum_N1e
    p2p <- rowSums(N2e * F_alt_Z) / sum_N2e

    # Pooled ALT variance — [2,2] element of the 2×2 covariance: n_s
    cm1_22 <- pm_alt * (1 - pm_alt) / (2 * N1_mat)
    cm2_22 <- pm_alt * (1 - pm_alt) / (2 * N2_mat)
    c1p_22 <- rowSums((cm1_22 + var_alt_C) * N1e^2) / sum_N1e^2
    c2p_22 <- rowSums((cm2_22 + var_alt_Z) * N2e^2) / sum_N2e^2

    # Wald test df=1: (delta_f)^2 / combined_var — all SNPs at once
    tstat       <- (p1p - p2p)^2 / (c1p_22 + c2p_22)
    Wald_log10p <- -log10(pchisq(tstat, 1L, lower.tail = FALSE))

    # ── Per-sample imputed ALT frequencies (long format) ─────────────────────
    means_tib <- bind_rows(
      tibble(chr = mychr, pos = rep(snps_i$POS, each = nrepl),
             TRT = "C", REP = rep(seq_len(nrepl), n_s),
             F_alt = as.vector(t(F_alt_C)), cM = rep(snps_i$cM, each = nrepl)),
      tibble(chr = mychr, pos = rep(snps_i$POS, each = nrepl),
             TRT = "Z", REP = rep(seq_len(nrepl), n_s),
             F_alt = as.vector(t(F_alt_Z)), cM = rep(snps_i$cM, each = nrepl))
    )

    list(
      scan = tibble(chr                    = mychr,
                    pos                    = snps_i$POS,
                    Wald_log10p            = Wald_log10p,
                    cM                     = snps_i$cM,
                    n_informative_founders = as.integer(rowSums(S > 0.5))),
      means = means_tib
    )
  })

results <- map(all_results, "scan")  %>% compact() %>% bind_rows() %>% filter(!is.na(Wald_log10p))
means   <- map(all_results, "means") %>% compact() %>% bind_rows()

cat("Writing", fileout, "\n")
write.table(results, fileout)
cat(sprintf("  %d SNPs\n", nrow(results)))

cat("Writing", fileout_means, "\n")
write.table(means, fileout_means)
cat(sprintf("  %d rows\n", nrow(means)))

cat("Done.\n")
