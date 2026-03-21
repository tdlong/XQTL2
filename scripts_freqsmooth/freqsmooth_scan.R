#!/usr/bin/env Rscript
###############################################################################
# freqsmooth_scan.R  —  Step 2 of the freqsmooth pipeline
#
# Reads the smoothed RDS from smooth_haps.R. For each window, reconstructs the
# per-replicate frequency matrices and covariance arrays from the long-format
# smoothed data, pools replicates with effective-N weighting, and runs the
# stabilized Wald test (eigenvalue-regularized). Heritability is computed from
# the same smoothed frequencies as the Wald test.
#
# Usage:
#   Rscript scripts_freqsmooth/freqsmooth_scan.R \
#       --chr chrX --dir process/proj/SCAN_NAME \
#       --outdir SCAN_NAME --rfile helpfiles/proj/design.txt
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
    "--chr"    = { parsed$chr    <- args[i+1]; i <- i+2L },
    "--dir"    = { parsed$dir    <- args[i+1]; i <- i+2L },
    "--outdir" = { parsed$outdir <- args[i+1]; i <- i+2L },
    "--rfile"  = { parsed$rfile  <- args[i+1]; i <- i+2L },
    stop(paste("Unknown argument:", args[i]))
  )
}

mychr <- parsed$chr
design.df <- read.table(parsed$rfile, header = TRUE)
source("scripts/scan_functions.R")

filein  <- file.path(parsed$dir, paste0(parsed$outdir, ".smooth.", mychr, ".rds"))
fileout <- file.path(parsed$dir, paste0(parsed$outdir, ".scan.", mychr, ".txt"))

options(dplyr.summarise.inform = FALSE)

# ── Load smoothed data ────────────────────────────────────────────────────────
cat("Reading", filein, "\n")
sm            <- readRDS(filein)
freq_smoothed <- sm$freq       # (CHROM, pos, TRT, REP, founder, freq, Num)
err_smoothed  <- sm$err        # (CHROM, pos, TRT, REP, fi, fj, v)
founder_names <- sm$founder_names
nrepl         <- sm$nrepl
nF            <- length(founder_names)

ProportionSelect <- design.df %>%
  filter(TRT == "Z") %>% select(REP, Proportion) %>% arrange(REP)

# ── Per-window scan using group_map ──────────────────────────────────────────
# group_map applies a function to each (CHROM, pos) group, returning a list
# that we bind into the final result table.
cat(sprintf("Running Wald scan on %s...\n", mychr))

scan_results <- freq_smoothed %>%
  group_by(CHROM, pos) %>%
  group_map(function(w_freq, key) {

    # ── Reconstruct p1, p2 matrices (nrepl x nF) ─────────────────────────────
    to_mat <- function(trt) {
      w_freq %>%
        filter(TRT == trt) %>%
        arrange(REP) %>%
        pivot_wider(names_from = founder, values_from = freq) %>%
        select(all_of(founder_names)) %>%
        as.matrix()
    }
    p1 <- to_mat("C"); p2 <- to_mat("Z")
    if (anyNA(p1) || anyNA(p2)) return(NULL)

    N1 <- w_freq %>% filter(TRT == "C") %>% arrange(REP) %>%
            group_by(REP) %>% slice(1) %>% pull(Num)
    N2 <- w_freq %>% filter(TRT == "Z") %>% arrange(REP) %>%
            group_by(REP) %>% slice(1) %>% pull(Num)

    # ── Reconstruct covar arrays (nF x nF x nrepl) ───────────────────────────
    w_err <- err_smoothed %>%
      filter(CHROM == key$CHROM, pos == key$pos)

    to_arr <- function(trt) {
      arr <- array(NA_real_, c(nF, nF, nrepl))
      for (r in seq_len(nrepl)) {
        mat_vals <- w_err %>% filter(TRT == trt, REP == r) %>%
                    arrange(fi, fj) %>% pull(v)
        if (length(mat_vals) == nF^2) arr[,,r] <- matrix(mat_vals, nF, nF)
      }
      arr
    }
    covar1 <- to_arr("C"); covar2 <- to_arr("Z")
    if (anyNA(covar1) || anyNA(covar2)) return(NULL)

    # ── Pool replicates (effective-N weighting) ───────────────────────────────
    N1e <- N2e <- numeric(nrepl)
    cv1 <- cv2 <- array(0, c(nF, nF, nrepl))
    for (r in seq_len(nrepl)) {
      pm  <- (N1[r]*p1[r,] + N2[r]*p2[r,]) / (N1[r]+N2[r])
      cm1 <- mn.covmat(pm, 2*N1[r]); cm2 <- mn.covmat(pm, 2*N2[r])
      N1e[r] <- sum(diag(cm1))*4*N1[r]^2 /
        (sum(diag(cm1))*2*N1[r] + 2*N1[r]*sum(diag(covar1[,,r])))
      N2e[r] <- sum(diag(cm2))*4*N2[r]^2 /
        (sum(diag(cm2))*2*N2[r] + 2*N2[r]*sum(diag(covar2[,,r])))
      cv1[,,r] <- (cm1 + covar1[,,r]) * N1e[r]^2
      cv2[,,r] <- (cm2 + covar2[,,r]) * N2e[r]^2
    }
    p1p <- as.vector(N1e %*% p1 / sum(N1e))
    p2p <- as.vector(N2e %*% p2 / sum(N2e))
    c1p <- rowSums(cv1, dims=2) / sum(N1e)^2
    c2p <- rowSums(cv2, dims=2) / sum(N2e)^2

    # ── Stabilized Wald test ──────────────────────────────────────────────────
    df    <- nF - 1
    covar <- c1p + c2p
    eg    <- eigen(covar)
    eval  <- pmax(eg$values[1:df], RIDGE_FRACTION * mean(eg$values[1:df]))
    trafo <- diag(1/sqrt(eval)) %*% t(eg$vectors[,1:df])
    tstat <- sum((trafo %*% (p1p - p2p))^2)
    pval  <- exp(pchisq(tstat, df, lower.tail=FALSE, log.p=TRUE))

    # ── Heritability (same smoothed p1, p2 as Wald test) ─────────────────────
    h2 <- tryCatch(
      Heritability(p1, p2, nrepl, ProportionSelect, AF_CUTOFF),
      error = function(e) list(Falconer_H2 = NA_real_, Cutler_H2 = NA_real_)
    )

    tibble(chr         = key$CHROM,
           pos         = key$pos,
           Wald_log10p = -log10(exp(pval)),
           Falc_H2     = h2$Falconer_H2,
           Cutl_H2     = h2$Cutler_H2)
  }) %>%
  bind_rows()

# ── Add genetic position and write ───────────────────────────────────────────
cat("Writing", fileout, "\n")
scan_results %>% add_genetic() %>% write.table(fileout)
cat("Done.\n")
