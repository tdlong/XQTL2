#!/usr/bin/env Rscript
###############################################################################
# plot_pseudoscan.R — 5-panel Wald -log10(p) line plot from haplotype scan
#
# Usage:
#   Rscript scripts/plot_pseudoscan.R \
#       --scan  process/proj/SCAN/SCAN.scan.txt \
#       --out   process/proj/SCAN/wald.png \
#       --format powerpoint
#
# Multi-scan overlay:
#   Rscript scripts/plot_pseudoscan.R \
#       --scan process/proj/SCAN_M/SCAN_M.scan.txt \
#       --scan process/proj/SCAN_F/SCAN_F.scan.txt \
#       --label Male --label Female \
#       --colour "#1F78B4" --colour "#E31A1C" \
#       --out overlay.png --format powerpoint
#
# Optional:
#   --threshold 10        dashed horizontal line
#   --genes genes.txt     tab-delimited: name, chr, pos_mb
#   --peaks peaks.txt     tab-delimited: label, chr, pos_mb
#   --height 7.0          override height in inches (default: 1.4 per chr)
###############################################################################

library(tidyverse)

# ── Parse command-line arguments ─────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

scan_files <- character(0)
scan_labels <- character(0)
scan_colours <- character(0)
out_file <- NULL
format_name <- "powerpoint"
threshold <- NULL
genes_file <- NULL
peaks_file <- NULL
out_height <- NULL

i <- 1
while (i <= length(args)) {
  switch(args[i],
    "--scan"      = { scan_files   <- c(scan_files,   args[i+1]); i <- i + 2 },
    "--label"     = { scan_labels  <- c(scan_labels,  args[i+1]); i <- i + 2 },
    "--colour"    = { scan_colours <- c(scan_colours, args[i+1]); i <- i + 2 },
    "--color"     = { scan_colours <- c(scan_colours, args[i+1]); i <- i + 2 },
    "--out"       = { out_file     <- args[i+1]; i <- i + 2 },
    "--format"    = { format_name  <- args[i+1]; i <- i + 2 },
    "--threshold" = { threshold    <- as.numeric(args[i+1]); i <- i + 2 },
    "--genes"     = { genes_file   <- args[i+1]; i <- i + 2 },
    "--peaks"     = { peaks_file   <- args[i+1]; i <- i + 2 },
    "--height"    = { out_height   <- as.numeric(args[i+1]); i <- i + 2 },
    stop("Unknown argument: ", args[i])
  )
}

if (length(scan_files) == 0) stop("At least one --scan is required")
if (is.null(out_file))       stop("--out is required")

if (length(scan_labels) == 0)  scan_labels  <- basename(scan_files)
if (length(scan_colours) == 0) scan_colours <- if (length(scan_files) == 1) "#1F78B4" else hcl.colors(length(scan_files), "Dark 2")

# ── Format presets ───────────────────────────────────────────────────────────
FORMAT_PRESETS <- list(
  manuscript_half       = c(w = 3.5, dpi = 300, font =  7),
  manuscript_full       = c(w = 7.0, dpi = 300, font =  8),
  manuscript_half_hires = c(w = 3.5, dpi = 600, font =  7),
  manuscript_full_hires = c(w = 7.0, dpi = 600, font =  8),
  powerpoint            = c(w = 8.0, dpi = 150, font =  9),
  web                   = c(w = 7.0, dpi = 150, font = 10),
  email                 = c(w = 6.0, dpi = 100, font =  9)
)

if (!format_name %in% names(FORMAT_PRESETS))
  stop("Unknown format '", format_name, "'. Choose from: ",
       paste(names(FORMAT_PRESETS), collapse = ", "))

fmt       <- FORMAT_PRESETS[[format_name]]
BASE_FONT <- fmt["font"]

# ── Constants ────────────────────────────────────────────────────────────────
chr_order  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
chr_labels <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")

# dm6 euchromatin boundaries
# Source: Huynh et al. 2023 PLoS Genet 19:e1010439, Supplementary Table S2
HET_BOUNDS <- read.table("helpfiles/het_bounds.txt", header = TRUE, comment.char = "#")

# ── Read scans ───────────────────────────────────────────────────────────────
scans_df <- map_dfr(seq_along(scan_files), function(j) {
  f <- scan_files[j]
  if (!file.exists(f)) { warning("File not found, skipping: ", f); return(NULL) }
  read.table(f, header = TRUE) %>%
    as_tibble() %>%
    mutate(pos_mb = pos / 1e6, label = scan_labels[j])
}) %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  filter(!is.na(Wald_log10p), chr %in% chr_order)

if (nrow(scans_df) == 0) stop("No scan data loaded — check --scan paths.")

if (is.null(out_height))
  out_height <- length(unique(as.character(scans_df$chr))) * 1.4

colour_map <- setNames(scan_colours, scan_labels)

chr_label_df <- scans_df %>%
  distinct(chr) %>%
  mutate(label = chr_labels[as.character(chr)])

# ── Heterochromatin shading ──────────────────────────────────────────────────
het_rects <- HET_BOUNDS %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  left_join(
    scans_df %>% group_by(chr) %>% summarise(xmax_data = max(pos_mb, na.rm = TRUE)),
    by = "chr"
  ) %>%
  { bind_rows(
      transmute(., chr, xmin = 0,       xmax = eu_start),
      transmute(., chr, xmin = eu_end,  xmax = xmax_data)
  )}

# ── Plot ─────────────────────────────────────────────────────────────────────
p <- ggplot(scans_df,
            aes(x = pos_mb, y = Wald_log10p, colour = label, group = label)) +
  geom_rect(data = het_rects,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.5, colour = NA, inherit.aes = FALSE) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ chr, ncol = 1, scales = "free") +
  geom_text(data = chr_label_df,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5,
            size = BASE_FONT * 0.25, colour = "grey30", fontface = "bold",
            inherit.aes = FALSE) +
  scale_colour_manual(values = colour_map, name = NULL) +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Position (Mb)", y = expression(-log[10](italic(P)) ~ "Wald")) +
  theme_classic(base_size = BASE_FONT) +
  theme(
    legend.position        = if (length(scan_files) > 1) "top" else "none",
    strip.background       = element_blank(),
    strip.text             = element_blank(),
    panel.grid.major.y     = element_line(colour = "grey92", linewidth = 0.3),
    panel.spacing          = unit(0.4, "lines")
  )

if (!is.null(threshold))
  p <- p + geom_hline(yintercept = threshold, linetype = "dashed",
                      colour = "grey40", linewidth = 0.4)

# ── Gene labels ──────────────────────────────────────────────────────────────
if (!is.null(genes_file) && file.exists(genes_file)) {
  genes_df <- read.delim(genes_file, stringsAsFactors = FALSE) %>%
    mutate(chr = factor(chr, levels = chr_order)) %>%
    filter(chr %in% chr_order)

  if (nrow(genes_df) > 0) {
    p <- p +
      geom_vline(data = genes_df, aes(xintercept = pos_mb),
                 colour = "darkgreen", linetype = "dotted",
                 linewidth = 0.5, inherit.aes = FALSE) +
      geom_text(data = genes_df, aes(x = pos_mb, y = Inf, label = name),
                angle = 90, hjust = 1.1, vjust = -0.5,
                size = BASE_FONT * 0.30, colour = "darkgreen",
                fontface = "italic", inherit.aes = FALSE)
  }
}

# ── Peak labels ──────────────────────────────────────────────────────────────
if (!is.null(peaks_file) && file.exists(peaks_file)) {
  peaks_df <- read.delim(peaks_file, stringsAsFactors = FALSE) %>%
    mutate(chr = factor(chr, levels = chr_order)) %>%
    filter(chr %in% chr_order)

  if (nrow(peaks_df) > 0) {
    # Find the max Wald across ALL scans at each peak's nearest position.
    # Loop avoids rowwise+sapply closure issues with .data pronoun.
    y_vals <- numeric(nrow(peaks_df))
    for (k in seq_len(nrow(peaks_df))) {
      pk_chr <- peaks_df$chr[k]
      pk_pos <- peaks_df$pos_mb[k]
      s <- scans_df %>% filter(chr == pk_chr)
      if (nrow(s) == 0) { y_vals[k] <- 0; next }
      best <- 0
      for (lbl in unique(s$label)) {
        ss <- s %>% filter(label == lbl)
        val <- ss$Wald_log10p[which.min(abs(ss$pos_mb - pk_pos))]
        if (val > best) best <- val
      }
      y_vals[k] <- best
    }
    peaks_df$y_val <- y_vals

    # Per-chromosome y ranges so offsets scale to each panel's local axis
    chr_ymax <- scans_df %>%
      group_by(chr) %>%
      summarise(y_range = max(Wald_log10p, na.rm = TRUE), .groups = "drop")

    peaks_df <- peaks_df %>%
      left_join(chr_ymax, by = "chr") %>%
      mutate(y_triangle = y_val + y_range * 0.04,
             y_label    = y_val + y_range * 0.12)

    p <- p +
      geom_point(data = peaks_df,
                 aes(x = pos_mb, y = y_triangle),
                 shape = 25, size = 2, fill = NA,
                 colour = "black", stroke = 0.8, inherit.aes = FALSE) +
      geom_text(data = peaks_df,
                aes(x = pos_mb, y = y_label, label = label),
                size = BASE_FONT * 0.30, colour = "black",
                hjust = 0.5, vjust = 0, inherit.aes = FALSE)
  }
}

# ── Save ─────────────────────────────────────────────────────────────────────
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
png(out_file, width = fmt["w"], height = out_height,
    units = "in", res = fmt["dpi"])
print(p)
dev.off()
cat("Saved:", out_file, "\n")
