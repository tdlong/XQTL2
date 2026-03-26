#!/usr/bin/env Rscript
###############################################################################
# plot_single_scan.R
#
# 5-panel genome scan plot (single trace) on the cM scale.
#
# Usage:
#   Rscript plot_single_scan.R <scan_dir> <outfile_prefix> <flymap> [label] [gene_file]
#
# Arguments:
#   scan_dir        Directory containing per-chromosome scan files
#   outfile_prefix  Output prefix (no extension); .png written
#   flymap          Path to flymap.r6.txt (Mb->cM table)
#   label           (Optional) trace label for the plot title, e.g. "AGE_SY20_M"
#   gene_file       (Optional) TSV with columns: gene_name, chr, pos_bp
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript plot_single_scan.R <scan_dir> <outfile_prefix> <flymap> [label] [gene_file]\n")
  quit(status = 1)
}

scan_dir    <- args[1]
outprefix   <- args[2]
flymap_file <- args[3]
label       <- if (length(args) >= 4) args[4] else basename(scan_dir)
gene_file   <- if (length(args) >= 5) args[5] else NULL

# ── Column names — auto-detected from file; override here if needed ───────────
COL_POS  <- "pos"   # physical position column (confirmed in pseudoN scan output)
# Stat column: try candidates in order; first match wins
STAT_CANDIDATES <- list(
  list(col = "Wald_log10p", is_log10p = TRUE),   # ZINC2-style pseudoN scan
  list(col = "log10p",      is_log10p = TRUE),    # malathion-style pseudoN scan
  list(col = "wald_stat",   is_log10p = FALSE),   # freqsmooth scan (chi-sq)
  list(col = "lrt_stat",    is_log10p = FALSE)
)

CHROMS <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")

COL_TRACE <- "#1F78B4"   # blue
LWD       <- 1.5

N_TESTS    <- 1e6
SIG_THRESH <- -log10(0.05 / N_TESTS)

FLYMAP_CHR <- "chr"
FLYMAP_BP  <- "pos_bp"
FLYMAP_CM  <- "pos_cM"

# ── Locate scan file for one chromosome ──────────────────────────────────────
scan_file <- function(dir, chr) {
  # Try fixed names first, then glob for *chr*.txt / *chr*.tsv
  fixed <- c(
    file.path(dir, paste0(chr, ".scan.tsv")),
    file.path(dir, paste0(chr, ".tsv")),
    file.path(dir, paste0(chr, ".txt")),
    file.path(dir, paste0("scan_", chr, ".tsv")),
    file.path(dir, paste0(chr, "_scan.tsv"))
  )
  found <- fixed[file.exists(fixed)]
  if (length(found) > 0) return(found[1])
  # Glob fallback: any file containing the chr name, excluding RefAlt/means/badloci
  globs <- Sys.glob(file.path(dir, paste0("*", chr, "*.txt")))
  globs <- c(globs, Sys.glob(file.path(dir, paste0("*", chr, "*.tsv"))))
  globs <- globs[!grepl("RefAlt|means|badloci|R\\.haps", globs)]
  if (length(globs) > 0) return(globs[1])
  NULL
}

# ── Read flymap ───────────────────────────────────────────────────────────────
message("Reading flymap: ", flymap_file)
flymap <- read.table(flymap_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ── Read one scan file ────────────────────────────────────────────────────────
read_scan <- function(filepath, chr) {
  d <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!COL_POS %in% colnames(d))
    stop("Column '", COL_POS, "' not found in ", filepath,
         "\n  Available: ", paste(colnames(d), collapse = ", "))
  # Auto-detect stat column
  stat_match <- Filter(function(x) x$col %in% colnames(d), STAT_CANDIDATES)
  if (length(stat_match) == 0)
    stop("No recognised stat column in ", filepath,
         "\n  Available: ", paste(colnames(d), collapse = ", "))
  COL_STAT_USE  <- stat_match[[1]]$col
  IS_LOG10P_USE <- stat_match[[1]]$is_log10p
  message("    stat column: ", COL_STAT_USE, "  (already log10p: ", IS_LOG10P_USE, ")")
  if (IS_LOG10P_USE) {
    d$neglog10p <- d[[COL_STAT_USE]]
  } else {
    d$neglog10p <- -pchisq(d[[COL_STAT_USE]], df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  }
  # Use cM column if present, otherwise interpolate from flymap
  if ("cM" %in% colnames(d)) {
    d$pos_cM <- d$cM
  } else {
    fm <- flymap[flymap[[FLYMAP_CHR]] == chr, ]
    if (nrow(fm) == 0) stop("No flymap rows for: ", chr)
    d$pos_cM <- approx(fm[[FLYMAP_BP]], fm[[FLYMAP_CM]], xout = d[[COL_POS]], rule = 2)$y
  }
  d
}

# ── Gene annotations ──────────────────────────────────────────────────────────
genes <- NULL
if (!is.null(gene_file) && file.exists(gene_file)) {
  genes <- read.table(gene_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  genes$pos_cM <- mapply(function(chr, bp) {
    fm <- flymap[flymap[[FLYMAP_CHR]] == chr, ]
    if (nrow(fm) == 0) return(NA_real_)
    approx(fm[[FLYMAP_BP]], fm[[FLYMAP_CM]], xout = bp, rule = 2)$y
  }, genes$chr, genes$pos_bp)
}

# ── Plot ──────────────────────────────────────────────────────────────────────
make_plot <- function() {
  par(mfrow = c(1, 5), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (chr in CHROMS) {
    message("  Plotting: ", chr)
    sf <- scan_file(scan_dir, chr)
    if (is.null(sf)) {
      warning("No scan file for ", chr, " — empty panel")
      plot(NA, xlim = c(0,1), ylim = c(0,1),
           xlab = "cM", ylab = expression(-log[10](p)), main = chr, las = 1, bty = "l")
      next
    }
    d <- read_scan(sf, chr)
    ylim_max <- max(c(d$neglog10p, SIG_THRESH), na.rm = TRUE) * 1.05
    plot(NA, xlim = range(d$pos_cM, na.rm = TRUE), ylim = c(0, ylim_max),
         xlab = "cM", ylab = expression(-log[10](p)), main = chr, las = 1, bty = "l")
    lines(d$pos_cM, d$neglog10p, col = COL_TRACE, lwd = LWD)
    abline(h = SIG_THRESH, lty = 2, col = "grey40", lwd = 1)
    if (!is.null(genes)) {
      cg <- genes[genes$chr == chr & !is.na(genes$pos_cM), ]
      for (i in seq_len(nrow(cg))) {
        abline(v = cg$pos_cM[i], lty = 2, lwd = 1.5, col = "#333333")
        text(cg$pos_cM[i], ylim_max * 0.97, cg$gene_name[i],
             srt = 90, adj = c(1, -0.3), cex = 0.9, col = "#333333")
      }
    }
  }
  mtext(paste("Genome scan —", label), outer = TRUE, cex = 1.2, font = 2)
}

png_file <- paste0(outprefix, ".png")
message("Writing PNG: ", png_file)
png(png_file, width = 18, height = 5, units = "in", res = 150)
make_plot()
invisible(dev.off())
message("Done: ", png_file)
