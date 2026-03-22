#!/usr/bin/env Rscript
###############################################################################
# plot_H2_overlay.R — 5-panel Falconer + Cutler H² overlay from haplotype scan
#
# Usage:
#   Rscript scripts/plot_H2_overlay.R \
#       --scan   process/proj/SCAN/SCAN.scan.txt \
#       --out    process/proj/SCAN/H2.png \
#       --format powerpoint
#
# Optional:
#   --threshold 0.5       dashed horizontal line
#   --genes genes.txt     tab-delimited: name, chr, pos_mb
#   --peaks peaks.txt     tab-delimited: label, chr, pos_mb
#   --height 7.0          override height in inches (default: 1.4 per chr)
#   --falc-colour "#E31A1C"   Falconer H² colour
#   --cutl-colour "#33A02C"   Cutler H² colour
###############################################################################

library(tidyverse)

# ── Parse command-line arguments ─────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

scan_file   <- NULL
out_file    <- NULL
format_name <- "powerpoint"
threshold   <- NULL
genes_file  <- NULL
peaks_file  <- NULL
out_height  <- NULL
falc_colour <- "#E31A1C"
cutl_colour <- "#33A02C"

i <- 1
while (i <= length(args)) {
  switch(args[i],
    "--scan"        = { scan_file   <- args[i+1]; i <- i + 2 },
    "--out"         = { out_file    <- args[i+1]; i <- i + 2 },
    "--format"      = { format_name <- args[i+1]; i <- i + 2 },
    "--threshold"   = { threshold   <- as.numeric(args[i+1]); i <- i + 2 },
    "--genes"       = { genes_file  <- args[i+1]; i <- i + 2 },
    "--peaks"       = { peaks_file  <- args[i+1]; i <- i + 2 },
    "--height"      = { out_height  <- as.numeric(args[i+1]); i <- i + 2 },
    "--falc-colour" = { falc_colour <- args[i+1]; i <- i + 2 },
    "--cutl-colour" = { cutl_colour <- args[i+1]; i <- i + 2 },
    stop("Unknown argument: ", args[i])
  )
}

if (is.null(scan_file)) stop("--scan is required")
if (is.null(out_file))  stop("--out is required")

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

HET_BOUNDS <- tribble(
  ~chr,    ~eu_start, ~eu_end,
  "chrX",   2.5,      21.2,
  "chr2L",  0.5,      22.9,
  "chr2R",  1.3,      25.1,
  "chr3L",  0.7,      24.0,
  "chr3R",  4.5,      32.0
)

# ── Read data ────────────────────────────────────────────────────────────────
if (!file.exists(scan_file)) stop("File not found: ", scan_file)

scan_df <- read.table(scan_file, header = TRUE) %>%
  as_tibble() %>%
  mutate(pos_mb = pos / 1e6,
         chr    = factor(chr, levels = chr_order)) %>%
  filter(chr %in% chr_order) %>%
  pivot_longer(cols = c(Falc_H2, Cutl_H2),
               names_to  = "estimator",
               values_to = "H2") %>%
  filter(!is.na(H2)) %>%
  mutate(estimator = factor(estimator,
                            levels = c("Falc_H2", "Cutl_H2"),
                            labels = c("Falconer H\u00b2", "Cutler H\u00b2")))

if (nrow(scan_df) == 0) stop("No H2 data loaded — check --scan path.")

if (is.null(out_height))
  out_height <- length(unique(as.character(scan_df$chr))) * 1.4

colour_map <- c("Falconer H\u00b2" = falc_colour, "Cutler H\u00b2" = cutl_colour)

chr_label_df <- scan_df %>%
  distinct(chr) %>%
  mutate(label_text = chr_labels[as.character(chr)])

# ── Heterochromatin shading ──────────────────────────────────────────────────
het_rects <- HET_BOUNDS %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  left_join(
    scan_df %>% group_by(chr) %>% summarise(xmax_data = max(pos_mb, na.rm = TRUE)),
    by = "chr"
  ) %>%
  { bind_rows(
      transmute(., chr, xmin = 0,       xmax = eu_start),
      transmute(., chr, xmin = eu_end,  xmax = xmax_data)
  )}

# ── Plot ─────────────────────────────────────────────────────────────────────
p <- ggplot(scan_df, aes(x = pos_mb, y = H2, colour = estimator)) +
  geom_rect(data = het_rects,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.5, colour = NA, inherit.aes = FALSE) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  facet_wrap(~ chr, ncol = 1, scales = "free") +
  geom_text(data = chr_label_df,
            aes(x = Inf, y = Inf, label = label_text),
            hjust = 1.1, vjust = 1.5,
            size = BASE_FONT * 0.25, colour = "grey30", fontface = "bold",
            inherit.aes = FALSE) +
  scale_colour_manual(values = colour_map, name = NULL) +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Position (Mb)", y = "H\u00b2") +
  theme_classic(base_size = BASE_FONT) +
  theme(
    legend.position        = "top",
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
    # Find H2 value at each peak for label placement
    wide_df <- read.table(scan_file, header = TRUE) %>%
      as_tibble() %>%
      mutate(pos_mb = pos / 1e6, chr = factor(chr, levels = chr_order))

    peaks_df <- peaks_df %>%
      rowwise() %>%
      mutate(y_val = {
        s <- wide_df %>% filter(chr == .data$chr)
        if (nrow(s) == 0) 0 else max(s$Falc_H2[which.min(abs(s$pos_mb - pos_mb))],
                                       s$Cutl_H2[which.min(abs(s$pos_mb - pos_mb))],
                                       na.rm = TRUE)
      }) %>%
      ungroup()

    y_range <- max(scan_df$H2, na.rm = TRUE)
    peaks_df <- peaks_df %>%
      mutate(y_label = y_val + y_range * 0.08)

    p <- p +
      geom_point(data = peaks_df,
                 aes(x = pos_mb, y = y_val),
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
