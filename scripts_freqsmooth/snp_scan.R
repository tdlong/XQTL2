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
# The stabilized Wald test (df=1) and biallelic heritability follow.
#
# Usage:
#   Rscript scripts_freqsmooth/snp_scan.R \
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

filein  <- file.path(parsed$dir, paste0(parsed$outdir, ".smooth.", mychr, ".rds"))
fileout <- file.path(parsed$dir, paste0(parsed$outdir, ".snp_scan.", mychr, ".txt"))

options(dplyr.summarise.inform = FALSE)

# ── Load smoothed haplotype data ──────────────────────────────────────────────
cat("Reading", filein, "\n")
sm            <- readRDS(filein)
freq_smoothed <- sm$freq
err_smoothed  <- sm$err
nrepl         <- sm$nrepl
nF            <- length(founders)

ProportionSelect <- design.df %>%
  filter(TRT == "Z") %>% select(REP, Proportion) %>% arrange(REP)

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
         Falc_H2=numeric(), Cutl_H2=numeric(), cM=numeric(),
         n_informative_founders=integer()) %>%
    write.table(fileout)
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

results <- snp_df %>%
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

    # ── Vectorized SNP frequency computation ─────────────────────────────────
    # F_alt_C[s, r] = S[s,] · H_C[r,]  →  S %*% t(H_C)  (n_snps x nrepl)
    F_alt_C <- S %*% t(H_C)   # n_snps x nrepl
    F_alt_Z <- S %*% t(H_Z)

    # ── One Wald test per SNP ─────────────────────────────────────────────────
    # (small nrepl so the pooling loop is fast)
    Wald_log10p <- Falc_H2 <- Cutl_H2 <- rep(NA_real_, n_s)

    for (si in seq_len(n_s)) {
      s_vec <- S[si, ]
      S2    <- rbind(1 - s_vec, s_vec)   # 2 x nF  (REF row, ALT row)

      p1_mat <- cbind(1 - F_alt_C[si,], F_alt_C[si,])   # nrepl x 2
      p2_mat <- cbind(1 - F_alt_Z[si,], F_alt_Z[si,])

      # Propagate covariance: Sigma_SNP = S2 %*% Sigma_hap %*% t(S2)
      cov1 <- array(NA_real_, c(2, 2, nrepl))
      cov2 <- array(NA_real_, c(2, 2, nrepl))
      for (r in seq_len(nrepl)) {
        cov1[,,r] <- S2 %*% E_C[,,r] %*% t(S2)
        cov2[,,r] <- S2 %*% E_Z[,,r] %*% t(S2)
      }

      # Pool replicates
      N1e <- N2e <- numeric(nrepl)
      cv1 <- cv2 <- array(0, c(2, 2, nrepl))
      for (r in seq_len(nrepl)) {
        pm  <- (N1[r]*p1_mat[r,] + N2[r]*p2_mat[r,]) / (N1[r]+N2[r])
        cm1 <- mn.covmat(pm, 2*N1[r]); cm2 <- mn.covmat(pm, 2*N2[r])
        N1e[r] <- sum(diag(cm1))*4*N1[r]^2 /
          (sum(diag(cm1))*2*N1[r] + 2*N1[r]*sum(diag(cov1[,,r])))
        N2e[r] <- sum(diag(cm2))*4*N2[r]^2 /
          (sum(diag(cm2))*2*N2[r] + 2*N2[r]*sum(diag(cov2[,,r])))
        cv1[,,r] <- (cm1 + cov1[,,r]) * N1e[r]^2
        cv2[,,r] <- (cm2 + cov2[,,r]) * N2e[r]^2
      }
      p1p <- as.vector(N1e %*% p1_mat / sum(N1e))
      p2p <- as.vector(N2e %*% p2_mat / sum(N2e))
      c1p <- rowSums(cv1, dims=2) / sum(N1e)^2
      c2p <- rowSums(cv2, dims=2) / sum(N2e)^2

      # Wald test df=1
      covar <- c1p + c2p
      eg    <- eigen(covar)
      eval  <- pmax(eg$values[1], RIDGE_FRACTION * mean(eg$values[1]))
      trafo <- matrix(1/sqrt(eval)) %*% t(eg$vectors[,1,drop=FALSE])
      tstat <- sum((trafo %*% (p1p - p2p))^2)
      Wald_log10p[si] <- -log10(exp(pchisq(tstat, 1, lower.tail=FALSE, log.p=TRUE)))

      # Biallelic heritability: pass ALT frequency column as a 1-founder matrix
      h2 <- tryCatch(
        Heritability(matrix(p1_mat[,2], nrepl, 1),
                     matrix(p2_mat[,2], nrepl, 1),
                     nrepl, ProportionSelect, AF_CUTOFF),
        error = function(e) list(Falconer_H2=NA_real_, Cutler_H2=NA_real_)
      )
      Falc_H2[si] <- h2$Falconer_H2
      Cutl_H2[si] <- h2$Cutler_H2
    }

    tibble(chr                    = mychr,
           pos                    = snps_i$POS,
           Wald_log10p            = Wald_log10p,
           Falc_H2                = Falc_H2,
           Cutl_H2                = Cutl_H2,
           cM                     = snps_i$cM,
           n_informative_founders = as.integer(rowSums(S > 0.5)))
  }) %>%
  bind_rows() %>%
  filter(!is.na(Wald_log10p))

cat("Writing", fileout, "\n")
write.table(results, fileout)
cat(sprintf("Done. %d SNPs written.\n", nrow(results)))
