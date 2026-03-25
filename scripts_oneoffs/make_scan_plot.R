#!/usr/bin/env Rscript
###############################################################################
# make_scan_plot.R
#
# Slide-ready 5-panel stacked Manhattan from pseudoscan.txt files.
# Uses XQTL_Manhattan_5panel() directly for single scans.
# Replicates the same style for the ZINC2 M+F overlay.
#
# Usage — single scan (aging, pupation, malathion):
#   Rscript scripts/make_scan_plot.R <scan.txt> <out.png> [width] [height]
#
# Usage — overlay male + female (ZINC2):
#   Rscript scripts/make_scan_plot.R <male.txt> <out.png> [width] [height] <female.txt>
#
# Usage — single scan + candidate gene lines (malathion):
#   Rscript scripts/make_scan_plot.R <scan.txt> <out.png> [width] [height] "" <genes.txt>
#
# genes.txt: tab-separated with header, columns: gene  chr  Mb
###############################################################################

suppressPackageStartupMessages(library(tidyverse))
source("scripts/XQTL_plotting_functions.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript make_scan_plot.R <scan1.txt> <out.png> [width] [height] [scan2.txt] [genes.txt]\n")
  quit(status = 1)
}

file1   <- args[1]
outfile <- args[2]
png_w   <- if (length(args) >= 3 && nchar(args[3]) > 0) as.numeric(args[3]) else 7.6
png_h   <- if (length(args) >= 4 && nchar(args[4]) > 0) as.numeric(args[4]) else 6.0
file2   <- if (length(args) >= 5 && nchar(args[5]) > 0) args[5] else NULL
genes_f <- if (length(args) >= 6 && nchar(args[6]) > 0) args[6] else NULL

chr_order <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

read_scan <- function(f) {
  d <- read.table(f, header = TRUE, sep = "", stringsAsFactors = FALSE)
  d$chr <- factor(d$chr, levels = chr_order)
  d$Mb  <- d$pos / 1e6
  d
}

# ── Build plot ────────────────────────────────────────────────────────────────
if (!is.null(file2)) {

  # ZINC2 overlay: replicate XQTL_Manhattan_5panel style with M (blue) + F (red)
  d1 <- read_scan(file1);  d1$sex <- "Male"
  d2 <- read_scan(file2);  d2$sex <- "Female"
  dat <- bind_rows(d1, d2)
  dat$sex <- factor(dat$sex, levels = c("Male", "Female"))

  label_data <- dat %>%
    group_by(chr) %>%
    summarise(x = max(Mb), .groups = "drop")

  p <- ggplot(dat, aes(x = Mb, y = Wald_log10p, color = sex)) +
    geom_point(size = 0.25) +
    facet_wrap(~ chr, ncol = 1, scales = "free") +
    geom_text(data = label_data,
              aes(x = Inf, y = Inf, label = chr),
              hjust = 1.1, vjust = 1.1, size = 3,
              color = "black", inherit.aes = FALSE) +
    scale_color_manual(values = c(Male = "#1F78B4", Female = "#E31A1C")) +
    scale_y_continuous(limits = function(y) c(0, max(10, ceiling(max(y))))) +
    labs(x = "Mb", y = "-log10(p-value)") +
    theme_bw() +
    theme(
      panel.spacing    = unit(0.1, "lines"),
      strip.background = element_blank(),
      strip.text       = element_blank(),
      legend.position  = "none"
    )

} else {

  # Single scan — call XQTL_Manhattan_5panel directly
  d1 <- read_scan(file1)
  p  <- XQTL_Manhattan_5panel(d1, cM = FALSE)

}

# ── Candidate gene vertical lines (malathion only) ────────────────────────────
if (!is.null(genes_f)) {
  genes <- read.table(genes_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  genes$chr <- factor(genes$chr, levels = chr_order)
  p <- p +
    geom_vline(data = genes, aes(xintercept = Mb),
               linetype = "dashed", color = "black", linewidth = 0.6,
               inherit.aes = FALSE) +
    geom_text(data = genes, aes(x = Mb, label = gene),
              y = Inf, vjust = 1.5, hjust = -0.15,
              size = 3.5, color = "black", inherit.aes = FALSE)
}

# ── Write PNG ─────────────────────────────────────────────────────────────────
options(bitmapType = "cairo")
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
png(outfile, width = png_w, height = png_h, units = "in", res = 200)
print(p)
dev.off()
message("Written: ", outfile)
