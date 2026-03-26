#!/usr/bin/env Rscript
###############################################################################
# plot_for_slides.R
#
# 5-panel horizontal genome scan for PowerPoint slides.
# Modeled on XQTL_Manhattan_5panel() but with ncol=5 (horizontal layout).
# Optionally overlays two scans (e.g. male + female) in different colors.
#
# Usage — single scan:
#   Rscript scripts/plot_for_slides.R <pseudoscan.txt> <out.png>
#
# Usage — overlay (male + female):
#   Rscript scripts/plot_for_slides.R <male.txt> <out.png> <female.txt>
#
# Examples:
#   Rscript scripts/plot_for_slides.R \
#       process/ZINC2/ZINC2_M/ZINC2_M.pseudoscan.txt figures/slides/ZINC2_M.png
#
#   Rscript scripts/plot_for_slides.R \
#       process/ZINC2/ZINC2_M/ZINC2_M.pseudoscan.txt figures/slides/ZINC2_overlay.png \
#       process/ZINC2/ZINC2_F/ZINC2_F.pseudoscan.txt
###############################################################################

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript plot_for_slides.R <scan1.txt> <out.png> [scan2.txt]\n")
  quit(status = 1)
}

file1   <- args[1]
outfile <- args[2]
file2   <- if (length(args) >= 3) args[3] else NULL

chr_order <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

read_scan <- function(f) {
  message("Reading: ", f)
  df <- as_tibble(read.table(f))
  df$chr <- factor(df$chr, levels = chr_order)
  df$Mb  <- df$pos / 1e6
  df
}

# ── Single scan ───────────────────────────────────────────────────────────────
if (is.null(file2)) {

  df <- read_scan(file1)

  p <- ggplot(df, aes(x = cM, y = Wald_log10p)) +
    geom_point(color = "#003262", size = 0.25) +
    facet_wrap(~ chr, ncol = 1, scales = "free") +
    scale_y_continuous(limits = function(y) c(0, max(10, ceiling(max(y, na.rm = TRUE))))) +
    labs(x = "cM", y = expression(-log[10](italic(p)))) +
    theme_bw() +
    theme(
      panel.spacing      = unit(0.3, "lines"),
      strip.background   = element_blank(),
      strip.text         = element_text(size = 11, face = "bold"),
      axis.text          = element_text(size = 9),
      axis.title         = element_text(size = 11)
    )

# ── Overlay: two scans colored by group ──────────────────────────────────────
} else {

  df1 <- read_scan(file1)
  df2 <- read_scan(file2)
  df1$group <- "Male"
  df2$group <- "Female"
  df  <- bind_rows(df1, df2)
  df$group <- factor(df$group, levels = c("Male", "Female"))

  p <- ggplot(df, aes(x = cM, y = Wald_log10p, color = group)) +
    geom_point(size = 0.25, alpha = 0.8) +
    facet_wrap(~ chr, ncol = 1, scales = "free") +
    scale_color_manual(values = c(Male = "#1F78B4", Female = "#E31A1C"),
                       name = NULL) +
    scale_y_continuous(limits = function(y) c(0, max(10, ceiling(max(y, na.rm = TRUE))))) +
    labs(x = "cM", y = expression(-log[10](italic(p)))) +
    theme_bw() +
    theme(
      panel.spacing      = unit(0.3, "lines"),
      strip.background   = element_blank(),
      strip.text         = element_text(size = 11, face = "bold"),
      axis.text          = element_text(size = 9),
      axis.title         = element_text(size = 11),
      legend.position    = "top",
      legend.text        = element_text(size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
options(bitmapType = "cairo")
message("Writing: ", outfile)
png(outfile, width = 4, height = 7, units = "in", res = 150)
print(p)
dev.off()
message("Done.")
