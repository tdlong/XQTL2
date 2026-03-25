#!/usr/bin/env Rscript
###############################################################################
# plot_zinc_overlay.R
#
# Overlaid 5-panel genome scan plot: ZINC2 males vs females on the cM scale.
# Optionally annotates candidate genes with vertical dashed lines.
#
# Usage (positional arguments):
#   Rscript plot_zinc_overlay.R <male_dir> <female_dir> <outfile_prefix> <flymap> [gene_file]
#
# Arguments:
#   male_dir        Directory containing per-chromosome scan files for males
#                   e.g., process/ZINC2_v2/ZINC2_M_freqs250
#   female_dir      Directory containing per-chromosome scan files for females
#                   e.g., process/ZINC2_v2/ZINC2_F_freqs250
#   outfile_prefix  Output file prefix (without extension); .pdf and .png written
#                   e.g., analysis/figures/zinc_overlay
#   flymap          Path to Mb -> cM conversion table
#                   e.g., helpfiles/flymap.r6.txt
#   gene_file       (Optional) TSV with columns: gene_name, chr, pos_bp
#                   If provided, draws a vertical dashed line + gene label on
#                   the appropriate chromosome panel for each gene.
#                   e.g., helpfiles/malathion_genes.txt
#
# Output:
#   <outfile_prefix>.pdf  — 18 x 5 inch, 5-panel horizontal layout
#   <outfile_prefix>.png  — same layout at 300 dpi
#
# TODO: Verify column names in actual scan output files before running.
#       Run head -1 on one scan file and adjust COL_POS / COL_STAT below.
# TODO: Verify column names in flymap file (assumed: chr, pos_bp, pos_cM).
###############################################################################

# ── Malathion candidate genes (FlyBase r6 positions) ─────────────────────────
# These are provided here for reference; to use them pass a gene_file TSV.
#
# Malathion candidates (FlyBase r6 positions):
# Ace:    chr3R, 9,069,721 bp
# Cyp6g1: chr2R, 9,064,616 bp
# Mdr65:  chr2R, 9,163,757 bp
#
# To create gene_file:
# echo -e "gene_name\tchr\tpos_bp\nAce\tchr3R\t9069721\nCyp6g1\tchr2R\t9064616\nMdr65\tchr2R\t9163757" > helpfiles/malathion_genes.txt

# ── Parse command-line arguments ─────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 5) {
  cat("Usage: Rscript plot_zinc_overlay.R <male_dir> <female_dir> <outfile_prefix> <flymap> [gene_file]\n")
  cat("\n")
  cat("  male_dir        directory with per-chromosome male scan files\n")
  cat("  female_dir      directory with per-chromosome female scan files\n")
  cat("  outfile_prefix  output path prefix (no extension); .pdf and .png are written\n")
  cat("  flymap          path to flymap.r6.txt (Mb->cM table)\n")
  cat("  gene_file       (optional) TSV with columns: gene_name, chr, pos_bp\n")
  cat("                  draws a vertical dashed line + label on each chromosome panel\n")
  quit(status = 1)
}

male_dir    <- args[1]
female_dir  <- args[2]
outprefix   <- args[3]
flymap_file <- args[4]
gene_file   <- if (length(args) == 5) args[5] else NULL

# ── TODO: Verify these column names against actual scan file headers ──────────
# Run: head -1 <male_dir>/<some_scan_file>
# Expected columns from concat_Chromosome_Scans.Andreas.sh output:
COL_POS           <- "pos"           # confirmed from ZINC2 pseudoscan header
COL_STAT          <- "Wald_log10p"  # confirmed from ZINC2 pseudoscan header
STAT_IS_NEGLOG10P <- TRUE           # Wald_log10p is already -log10(p)

# ── Chromosomes to plot (in panel order) ─────────────────────────────────────
CHROMS <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")
# TODO: confirm chromosome name format matches scan file naming convention.
#       Files may be named e.g. scan_chr2L.tsv, 2L.scan.txt, etc.
#       Adjust SCAN_PATTERN below accordingly.

# ── TODO: Verify scan file naming convention ──────────────────────────────────
# Pattern used to locate per-chromosome scan files inside male_dir / female_dir.
# Adjust the glob pattern to match actual filenames.
# Examples:
#   paste0(chr, ".scan.tsv")
#   paste0("scan.", chr, ".txt")
#   paste0(chr, "_scan.tsv")
scan_file <- function(dir, chr) {
  # Try fixed naming patterns first
  candidates <- c(
    file.path(dir, paste0(chr, ".scan.tsv")),
    file.path(dir, paste0(chr, ".tsv")),
    file.path(dir, paste0(chr, ".txt")),
    file.path(dir, paste0("scan_", chr, ".tsv")),
    file.path(dir, paste0(chr, "_scan.tsv"))
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) return(found[1])
  # Glob fallback: any file containing the chr name (excludes RefAlt/means/haps)
  globs <- c(Sys.glob(file.path(dir, paste0("*", chr, "*.txt"))),
             Sys.glob(file.path(dir, paste0("*", chr, "*.tsv"))))
  globs <- globs[!grepl("RefAlt|means|badloci|R\\.haps", globs)]
  if (length(globs) > 0) return(globs[1])
  NULL
}

# ── Read flymap (Mb -> cM conversion) ────────────────────────────────────────
# TODO: Verify column names in flymap.r6.txt.
#       Assumed: chr (chromosome), pos_bp (physical position, bp),
#                pos_cM (genetic position, cM).
#       Adjust FLYMAP_CHR / FLYMAP_BP / FLYMAP_CM if different.
FLYMAP_CHR <- "chr"     # TODO: verify
FLYMAP_BP  <- "pos_bp"  # TODO: verify — may be "bp", "position_bp", "pos", etc.
FLYMAP_CM  <- "pos_cM"  # TODO: verify — may be "cM", "genetic_pos", etc.

message("Reading flymap: ", flymap_file)
flymap <- read.table(flymap_file, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
# Ensure flymap chromosome names match CHROMS vector above
# TODO: if flymap uses "2L" instead of "chr2L", prefix with "chr":
# flymap[[FLYMAP_CHR]] <- paste0("chr", flymap[[FLYMAP_CHR]])

# ── Helper: read one scan file and join with flymap ───────────────────────────
read_scan <- function(filepath, chr) {
  d <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Check expected columns exist
  if (!COL_POS %in% colnames(d)) {
    stop(paste0("Column '", COL_POS, "' not found in ", filepath,
                "\n  Available columns: ", paste(colnames(d), collapse = ", "),
                "\n  TODO: Update COL_POS at the top of this script."))
  }
  if (!COL_STAT %in% colnames(d)) {
    stop(paste0("Column '", COL_STAT, "' not found in ", filepath,
                "\n  Available columns: ", paste(colnames(d), collapse = ", "),
                "\n  TODO: Update COL_STAT at the top of this script."))
  }

  # Compute -log10(p) from Wald statistic (chi-sq df=1) if needed
  if (STAT_IS_NEGLOG10P) {
    d$neglog10p <- d[[COL_STAT]]
  } else {
    # Wald statistic assumed chi-sq distributed with 1 df
    # TODO: confirm degrees of freedom — adjust if multi-allelic (df > 1)
    d$neglog10p <- -pchisq(d[[COL_STAT]], df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  }

  # Join flymap for this chromosome to get cM positions
  fm <- flymap[flymap[[FLYMAP_CHR]] == chr, ]
  if (nrow(fm) == 0) {
    stop(paste0("No flymap rows found for chromosome: ", chr,
                "\n  TODO: Check chromosome name format in flymap vs CHROMS vector."))
  }

  # Interpolate cM from physical position using flymap
  d$pos_cM <- approx(x = fm[[FLYMAP_BP]], y = fm[[FLYMAP_CM]],
                     xout = d[[COL_POS]], rule = 2)$y
  d
}

# ── Read gene annotation file (optional) ─────────────────────────────────────
# Expected columns: gene_name, chr, pos_bp
# cM positions are interpolated from the flymap for each gene.
genes <- NULL
if (!is.null(gene_file)) {
  message("Reading gene file: ", gene_file)
  genes <- read.table(gene_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
  required_cols <- c("gene_name", "chr", "pos_bp")
  missing_cols  <- setdiff(required_cols, colnames(genes))
  if (length(missing_cols) > 0) {
    stop(paste0("gene_file is missing required column(s): ",
                paste(missing_cols, collapse = ", "),
                "\n  Expected columns: gene_name, chr, pos_bp"))
  }
  # Interpolate cM position for each gene using flymap
  genes$pos_cM <- mapply(function(chr, bp) {
    fm <- flymap[flymap[[FLYMAP_CHR]] == chr, ]
    if (nrow(fm) == 0) {
      warning(paste0("No flymap rows for gene chromosome: ", chr, " — skipping."))
      return(NA_real_)
    }
    approx(x = fm[[FLYMAP_BP]], y = fm[[FLYMAP_CM]], xout = bp, rule = 2)$y
  }, genes$chr, genes$pos_bp)
  message("  Loaded ", nrow(genes), " gene annotation(s).")
}

# ── Significance threshold ────────────────────────────────────────────────────
# Bonferroni threshold: -log10(0.05 / number_of_tests)
# TODO: replace with permutation-derived threshold once available.
N_TESTS    <- 1e6     # TODO: set to actual number of tested positions
SIG_THRESH <- -log10(0.05 / N_TESTS)

# ── Colors ────────────────────────────────────────────────────────────────────
COL_MALE   <- "#1F78B4"
COL_FEMALE <- "#E31A1C"
LWD        <- 1.5

# Gene annotation style
COL_GENE <- "#333333"  # dark gray
LTY_GENE <- 2          # dashed
LWD_GENE <- 1.5

# ── Plotting function ─────────────────────────────────────────────────────────
make_plot <- function() {
  par(mfrow = c(1, 5),
      mar   = c(4, 4, 3, 1),
      oma   = c(0, 0, 2, 0))

  for (chr in CHROMS) {
    message("  Plotting: ", chr)

    mfile <- scan_file(male_dir,   chr)
    ffile <- scan_file(female_dir, chr)

    # -- Safety: if scan file missing for this chromosome, draw empty panel --
    if (is.null(mfile) || is.null(ffile)) {
      warning(paste0("Scan file(s) not found for ", chr, " — plotting empty panel."))
      plot(NA,
           xlim = c(0, 1), ylim = c(0, 1),
           xlab = "cM", ylab = expression(-log[10](p)),
           main = chr,
           las  = 1, bty = "l")
      next
    }

    md <- read_scan(mfile, chr)
    fd <- read_scan(ffile, chr)

    ylim_max <- max(c(md$neglog10p, fd$neglog10p, SIG_THRESH), na.rm = TRUE) * 1.05
    xlim     <- range(c(md$pos_cM, fd$pos_cM), na.rm = TRUE)

    plot(NA,
         xlim = xlim, ylim = c(0, ylim_max),
         xlab = "cM", ylab = expression(-log[10](p)),
         main = chr,
         las  = 1, bty = "l")

    # Female trace (red, drawn first so male goes on top)
    lines(fd$pos_cM, fd$neglog10p, col = COL_FEMALE, lwd = LWD)

    # Male trace (blue)
    lines(md$pos_cM, md$neglog10p, col = COL_MALE, lwd = LWD)

    # Significance threshold
    abline(h = SIG_THRESH, lty = 2, col = "grey40", lwd = 1)

    # Gene annotations for this chromosome
    if (!is.null(genes)) {
      chr_genes <- genes[genes$chr == chr & !is.na(genes$pos_cM), ]
      if (nrow(chr_genes) > 0) {
        for (i in seq_len(nrow(chr_genes))) {
          gx   <- chr_genes$pos_cM[i]
          gnam <- chr_genes$gene_name[i]
          # Vertical dashed line
          abline(v = gx, lty = LTY_GENE, lwd = LWD_GENE, col = COL_GENE)
          # Gene label near the top of the panel, rotated 90 degrees
          text(x      = gx,
               y      = ylim_max * 0.97,
               labels = gnam,
               srt    = 90,
               adj    = c(1, -0.3),
               cex    = 0.9,
               col    = COL_GENE)
        }
      }
    }
  }

  # Overall title and legend (in outer margin)
  mtext("ZINC2 genome scan — males (blue) vs females (red)",
        outer = TRUE, cex = 1.2, font = 2)

  # Legend in first panel — drawn in the last call so it goes on top
  # Re-draw in the first panel by resetting to panel 1 is not straightforward
  # in base R; instead, add a legend to the figure margin.
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "")
  legend("topright",
         legend = c("Males", "Females"),
         col    = c(COL_MALE, COL_FEMALE),
         lwd    = LWD,
         bty    = "n",
         cex    = 0.9)
}

# ── Write PDF ─────────────────────────────────────────────────────────────────
pdf_file <- paste0(outprefix, ".pdf")
message("Writing PDF: ", pdf_file)
pdf(pdf_file, width = 18, height = 5)
make_plot()
invisible(dev.off())

# ── Write PNG ─────────────────────────────────────────────────────────────────
png_file <- paste0(outprefix, ".png")
message("Writing PNG: ", png_file)
png(png_file, width = 18, height = 5, units = "in", res = 300)
make_plot()
invisible(dev.off())

message("Done.")
message("  PDF: ", pdf_file)
message("  PNG: ", png_file)
