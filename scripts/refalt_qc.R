#!/usr/bin/env Rscript
###############################################################################
# refalt_qc.R — per-sample QC table from the RefAlt count tables
#
# One row per sample per chromosome. Every metric is computed from that sample's
# own counts; nothing is relative to other samples.
#
# Why this exists: a sample with poor coverage does not fail loudly in the
# current pipeline. est_hap2 fits each sample independently, so a bad BAM cannot
# corrupt its neighbours -- but between roughly 1% and 30% window coverage it
# returns haplotype frequencies that are badly wrong and still satisfy every
# constraint the pipeline checks. The Wald test partly protects itself (a large
# reconstruction covariance drives that replicate's effective N down), but
# heritability does not -- Falconer and Cutler H2 use the raw frequencies with no
# such weighting. Below ~1% coverage the estimate goes NA and hap_scan.R drops
# the window for every replicate. So: look before you scan.
#
# Columns
#   n_sites          catalog sites on this chromosome
#   median_depth     median REF+ALT. Median, not mean -- collapsed repeats
#                    inflate the mean and hide a thin library.
#   pct_zero         % of sites with no reads
#   pct_lt5          % of sites under 5 reads, where the per-SNP frequency is
#                    too noisy to inform the haplotype fit
#   pct_window_covered  covered sites as a % of what the catalog offers in that
#                    window. THE metric to read, and the one flagged on: it
#                    isolates this sample's deficit from how dense the catalog is.
#   sites_per_window median covered sites per haplotype window. THE absolute
#                    read: est_hap2 degrades as a function of covered sites per
#                    window, not of genome-wide depth, and clumped coverage can
#                    leave a sample underdetermined at respectable mean depth.
#                    Computed with the base half-window `size`; real windows
#                    expand in low-recombination regions, so this is the
#                    pessimistic case.
#   disp_k           robust overdispersion index. Negative-binomial size fitted
#                    by moments to depths within [0.25x, 4x] of the median.
#                    Poisson coverage would be very large; these libraries are
#                    never Poisson, but a bad one is far worse. Depth-invariant,
#                    unlike var/mean, which grows with mean depth and so cannot
#                    compare samples sequenced to different depths.
#                    RELATIVE, not calibrated: the band biases it upward for
#                    genuinely bad libraries (a true NB size of 1.0 reads ~2.3).
#                    Use it to rank and flag, not as an absolute NB size.
#   flag             OK / LOW_COVERAGE / PATCHY
#
# Usage:
#   Rscript scripts/refalt_qc.R \
#       --dir     process/<project> \
#       --parfile helpfiles/<project>/hap_params.R \
#       [--out    process/<project>/Calls/refalt_qc.txt]
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

# ── Flag thresholds ──────────────────────────────────────────────────────────
# sites_per_window is judged against the number of founders being fit. The rank
# condition (more equations than founders) is nowhere near sufficient: simulating
# est_hap2 on an 8-founder window shows the fit is already off by 0.15 in the
# worst founder at 73 covered sites and by 0.24 at 18, while still returning a
# result that satisfies every constraint the pipeline checks. Error falls below
# ~0.05 only in the high hundreds. So the floor is set well above the rank
# condition -- it is a "this fit is not trustworthy" line, not a "this fit is
# impossible" line.
# A sample is judged against what its own catalog makes available, not against an
# absolute site count: how densely the catalog covers the genome is a property of
# the project and affects every sample equally, so it is reported once rather than
# flagged per sample.
MIN_PCT_WINDOW_COVERED <- 50    # covered < half the window's sites -> LOW_COVERAGE
MAX_PCT_ZERO           <- 20    # more than this fraction uncovered -> LOW_COVERAGE
MIN_DISP_K             <- 3     # below this -> PATCHY
MIN_SITES_PER_FOUNDER  <- 25    # project-level: catalog too sparse for this nF

args   <- commandArgs(trailingOnly = TRUE)
parsed <- list()
i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--dir"     = { parsed$dir     <- args[i+1]; i <- i+2L },
    "--parfile" = { parsed$parfile <- args[i+1]; i <- i+2L },
    "--out"     = { parsed$out     <- args[i+1]; i <- i+2L },
    stop(paste("Unknown argument:", args[i]))
  )
}
for (need in c("dir","parfile"))
  if (is.null(parsed[[need]])) stop("Missing required argument: --", need)
if (is.null(parsed$out)) parsed$out <- file.path(parsed$dir, "Calls", "refalt_qc.txt")

source(parsed$parfile)
if (!exists("founders")) stop("No `founders` vector in ", parsed$parfile)
if (!exists("size") || !exists("step")) stop("No `size`/`step` in ", parsed$parfile)
nF <- length(founders)
cat(sprintf("%d founders, size=%d, step=%d\n", nF, size, step))

# ── Robust overdispersion ────────────────────────────────────────────────────
# Method of moments on a depth band. The band matters: a collapsed-repeat tail
# of 1% of sites at 20x depth drags an untrimmed estimate from a true 10 to 0.35,
# i.e. reports a good library as catastrophic.
disp_k <- function(d) {
  md <- median(d)
  if (!is.finite(md) || md <= 0) return(NA_real_)
  y <- d[d >= 0.25*md & d <= 4*md]
  if (length(y) < 100) return(NA_real_)
  m <- mean(y); v <- var(y)
  if (!is.finite(v) || v <= m) return(Inf)      # at or below Poisson
  m^2 / (v - m)
}

# Median covered sites per haplotype window, via a rolling count over the
# sorted covered positions.
sites_per_window <- function(pos_sel, all_pos) {
  if (length(pos_sel) == 0L) return(0)
  centres <- seq(min(all_pos), max(all_pos), by = step)
  if (length(centres) == 0L) return(0)
  p  <- sort(pos_sel)
  median(findInterval(centres + size, p) - findInterval(centres - size, p))
}

chrs <- c("chrX","chr2L","chr2R","chr3L","chr3R")
rows <- list()

for (mychr in chrs) {
  f <- file.path(parsed$dir, "Calls", paste0("RefAlt.", mychr, ".txt"))
  if (!file.exists(f)) { cat("skip (absent):", f, "\n"); next }
  cat("Reading", f, "...\n")
  ra <- fread(f, header = TRUE)

  samples <- unique(sub("^REF_", "", grep("^REF_", names(ra), value = TRUE)))
  # What the catalog makes available in a window, before any sample's coverage.
  avail <- sites_per_window(ra$POS, ra$POS)
  cat(sprintf("  %d sites, %d samples, %d catalog sites per window\n",
              nrow(ra), length(samples), avail))
  if (avail < MIN_SITES_PER_FOUNDER * nF)
    cat(sprintf("  NOTE: only %d catalog sites per %d bp window for %d founders.\n",
                avail, 2*size, nF),
        sprintf("        Every sample's haplotype fit is thinly determined here;\n"),
        sprintf("        this is the catalog, not any one sample.\n"), sep = "")

  for (s in samples) {
    rc <- ra[[paste0("REF_", s)]]; ac <- ra[[paste0("ALT_", s)]]
    if (is.null(rc) || is.null(ac)) next
    d <- as.numeric(rc) + as.numeric(ac)
    d[is.na(d)] <- 0

    spw <- sites_per_window(ra$POS[d > 0], ra$POS)
    pwc <- if (avail > 0) 100 * spw / avail else NA_real_
    k   <- disp_k(d)
    pz  <- 100 * mean(d == 0)

    flag <- "OK"
    if (is.finite(k) && k < MIN_DISP_K)                          flag <- "PATCHY"
    if ((is.finite(pwc) && pwc < MIN_PCT_WINDOW_COVERED) || pz > MAX_PCT_ZERO)
                                                                 flag <- "LOW_COVERAGE"

    rows[[length(rows)+1]] <- data.table(
      sample           = s,
      chr              = mychr,
      n_sites          = nrow(ra),
      median_depth     = median(d),
      pct_zero         = round(pz, 2),
      pct_lt5          = round(100 * mean(d < 5), 2),
      sites_per_window = spw,
      pct_window_covered = round(pwc, 1),
      disp_k           = round(k, 2),
      flag             = flag)
  }
}

if (length(rows) == 0L) stop("No RefAlt tables found under ", file.path(parsed$dir, "Calls"))

qc <- rbindlist(rows)
setorder(qc, sample, chr)
dir.create(dirname(parsed$out), showWarnings = FALSE, recursive = TRUE)
fwrite(qc, parsed$out, sep = "\t")
cat("\nWrote", parsed$out, "\n\n")

bad <- qc[flag != "OK"]
if (nrow(bad) == 0) {
  cat("All", uniqueN(qc$sample), "samples OK on all chromosomes.\n")
} else {
  cat(sprintf("%d sample x chromosome combination(s) flagged:\n\n", nrow(bad)))
  print(bad[, .(sample, chr, median_depth, pct_zero, sites_per_window,
                pct_window_covered, disp_k, flag)])
  cat("\nLOW_COVERAGE: the sample covers less than half the catalog sites its\n")
  cat("windows offer. Such a sample returns haplotype frequencies that can be\n")
  cat("badly wrong while still passing every check the pipeline makes, and it\n")
  cat("corrupts H2 even where the Wald test downweights it.\n")
  cat("PATCHY: coverage is far more clumped than its peers' -- reads are present\n")
  cat("but concentrated, so windows vary wildly in how well they are determined.\n")
}
