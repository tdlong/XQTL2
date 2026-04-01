#!/usr/bin/env Rscript
###############################################################################
# plot_manhattan.R — traditional single-row Manhattan plot from haplotype scan
#
# Chromosomes are concatenated on a single x-axis in fly order:
# chrX, chr2L, chr2R, chr3L, chr3R.  Heterochromatic regions are shaded
# continuously across chromosome boundaries.  Dotted vertical lines mark
# chromosome boundaries.
#
# Usage:
#   Rscript scripts/plot_manhattan.R \
#       --scan  process/proj/SCAN/SCAN.scan.txt \
#       --out   process/proj/SCAN/manhattan.png \
#       --format powerpoint
#
# Multi-scan overlay:
#   Rscript scripts/plot_manhattan.R \
#       --scan process/proj/SCAN_M/SCAN_M.scan.txt \
#       --scan process/proj/SCAN_F/SCAN_F.scan.txt \
#       --label Male --label Female \
#       --colour "#1F78B4" --colour "#E31A1C" \
#       --out overlay.png --format powerpoint
#
# Optional:
#   --column Cutl_H2      y-axis column (default: Wald_log10p)
#   --threshold 10        dashed horizontal line
#   --genes genes.txt     tab-delimited: name, chr, pos_mb
#   --peaks peaks.txt     tab-delimited: label, chr, pos_mb
#   --height 3.5          override height in inches (default: 3.5)
###############################################################################

library(tidyverse)

# ── Parse command-line arguments ─────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

scan_files <- character(0)
scan_labels <- character(0)
scan_colours <- character(0)
out_file <- NULL
format_name <- "powerpoint"
y_column <- "Wald_log10p"
threshold <- NULL
genes_file <- NULL
peaks_file <- NULL
out_height <- 3.5

i <- 1
while (i <= length(args)) {
  switch(args[i],
    "--scan"      = { scan_files   <- c(scan_files,   args[i+1]); i <- i + 2 },
    "--label"     = { scan_labels  <- c(scan_labels,  args[i+1]); i <- i + 2 },
    "--colour"    = { scan_colours <- c(scan_colours, args[i+1]); i <- i + 2 },
    "--color"     = { scan_colours <- c(scan_colours, args[i+1]); i <- i + 2 },
    "--out"       = { out_file     <- args[i+1]; i <- i + 2 },
    "--format"    = { format_name  <- args[i+1]; i <- i + 2 },
    "--column"    = { y_column     <- args[i+1]; i <- i + 2 },
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

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE))))
HET_BOUNDS <- read.table(file.path(script_dir, "../helpfiles/het_bounds.txt"), header = TRUE, comment.char = "#")

# ── Read scans ───────────────────────────────────────────────────────────────
scans_df <- map_dfr(seq_along(scan_files), function(j) {
  f <- scan_files[j]
  if (!file.exists(f)) { warning("File not found, skipping: ", f); return(NULL) }
  read.table(f, header = TRUE) %>%
    as_tibble() %>%
    mutate(pos_mb = pos / 1e6, label = scan_labels[j])
}) %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  filter(!is.na(.data[[y_column]]), chr %in% chr_order)

if (nrow(scans_df) == 0) stop("No scan data loaded — check --scan paths.")

# ── Build cumulative x-axis ─────────────────────────────────────────────────
chr_lengths <- scans_df %>%
  group_by(chr) %>%
  summarise(chr_max = max(pos_mb, na.rm = TRUE), .groups = "drop") %>%
  arrange(factor(chr, levels = chr_order))

chr_lengths$offset <- cumsum(c(0, head(chr_lengths$chr_max, -1)))
offsets <- setNames(chr_lengths$offset, chr_lengths$chr)

scans_df <- scans_df %>%
  mutate(x = pos_mb + offsets[as.character(chr)])

# Chromosome boundary positions (dotted vertical lines)
boundaries <- chr_lengths$offset[-1]

# Chromosome label positions (midpoint of each chromosome)
chr_mids <- chr_lengths %>%
  mutate(mid = offset + chr_max / 2,
         lbl = chr_labels[as.character(chr)])

# ── Heterochromatin shading (continuous across boundaries) ───────────────────
het_rects <- HET_BOUNDS %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  filter(chr %in% names(offsets)) %>%
  left_join(
    scans_df %>% group_by(chr) %>% summarise(xmax_data = max(pos_mb, na.rm = TRUE), .groups = "drop"),
    by = "chr"
  ) %>%
  mutate(off = offsets[as.character(chr)]) %>%
  { bind_rows(
      transmute(., xmin = off,           xmax = off + eu_start),
      transmute(., xmin = off + eu_end,  xmax = off + xmax_data)
  ) }

# ── Y-axis label ────────────────────────────────────────────────────────────
y_label <- if (y_column == "Wald_log10p") {
  expression(-log[10](italic(P)) ~ "Wald")
} else if (grepl("_H2$", y_column)) {
  bquote(italic(H)^2 ~ .(sub("_H2$", "", y_column)))
} else {
  y_column
}

colour_map <- setNames(scan_colours, scan_labels)

# ── Plot ─────────────────────────────────────────────────────────────────────
p <- ggplot(scans_df, aes(x = x, y = .data[[y_column]], colour = label, group = label)) +
  geom_rect(data = het_rects,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.5, colour = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = boundaries, linetype = "dotted", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = colour_map, name = NULL) +
  scale_x_continuous(
    breaks = chr_mids$mid,
    labels = chr_mids$lbl,
    expand = expansion(0)
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = y_label) +
  theme_classic(base_size = BASE_FONT) +
  theme(
    legend.position        = if (length(scan_files) > 1) "top" else "none",
    panel.grid.major.y     = element_line(colour = "grey92", linewidth = 0.3),
    axis.ticks.x           = element_blank()
  )

if (!is.null(threshold))
  p <- p + geom_hline(yintercept = threshold, linetype = "dashed",
                      colour = "grey40", linewidth = 0.4)

# ── Gene labels ──────────────────────────────────────────────────────────────
if (!is.null(genes_file) && file.exists(genes_file)) {
  genes_df <- read.delim(genes_file, stringsAsFactors = FALSE) %>%
    mutate(chr = factor(chr, levels = chr_order)) %>%
    filter(chr %in% chr_order) %>%
    mutate(x = pos_mb + offsets[as.character(chr)])

  if (nrow(genes_df) > 0) {
    p <- p +
      geom_vline(data = genes_df, aes(xintercept = x),
                 colour = "darkgreen", linetype = "dotted",
                 linewidth = 0.5, inherit.aes = FALSE) +
      geom_text(data = genes_df, aes(x = x, y = Inf, label = name),
                angle = 90, hjust = 1.1, vjust = -0.5,
                size = BASE_FONT * 0.30, colour = "darkgreen",
                fontface = "italic", inherit.aes = FALSE)
  }
}

# ── Peak labels ──────────────────────────────────────────────────────────────
if (!is.null(peaks_file) && file.exists(peaks_file)) {
  peaks_df <- read.delim(peaks_file, stringsAsFactors = FALSE) %>%
    mutate(chr = factor(chr, levels = chr_order)) %>%
    filter(chr %in% chr_order) %>%
    mutate(x = pos_mb + offsets[as.character(chr)])

  if (nrow(peaks_df) > 0) {
    y_vals <- numeric(nrow(peaks_df))
    for (k in seq_len(nrow(peaks_df))) {
      pk_x <- peaks_df$x[k]
      s <- scans_df %>% filter(abs(x - pk_x) == min(abs(x - pk_x)))
      y_vals[k] <- max(s[[y_column]], na.rm = TRUE)
    }
    peaks_df$y_val <- y_vals

    y_range <- max(scans_df[[y_column]], na.rm = TRUE)
    peaks_df <- peaks_df %>%
      mutate(y_triangle = y_val + y_range * 0.04,
             y_label    = y_val + y_range * 0.12)

    p <- p +
      geom_point(data = peaks_df,
                 aes(x = x, y = y_triangle),
                 shape = 25, size = 2, fill = NA,
                 colour = "black", stroke = 0.8, inherit.aes = FALSE) +
      geom_text(data = peaks_df,
                aes(x = x, y = y_label, label = label),
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
