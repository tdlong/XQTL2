#!/usr/bin/env Rscript
# MALATHION_TEST_v2_smooth125.R
#
# 4-panel Manhattan figure for the freqsmooth pipeline test:
#   Panel 1 — Haplotype Wald  -log10(p)
#   Panel 2 — SNP Wald  -log10(p)
#   Panel 3 — Falconer H²  (haplotype scan)
#   Panel 4 — Cutler H²   (haplotype scan)
#
# Run from XQTL2 project root:
#   Rscript helpfiles/malathion_test/MALATHION_TEST_v2_smooth125.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

SCAN_DIR <- "process/malathion_test/MALATHION_TEST_v2_smooth125"
SCAN     <- "MALATHION_TEST_v2_smooth125"
CHRS     <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
OUT_FILE <- file.path(SCAN_DIR, paste0(SCAN, ".4panel.png"))

# ── Load per-chromosome files ─────────────────────────────────────────────────
read_chr_files <- function(suffix) {
  map_dfr(CHRS, function(chr) {
    f <- file.path(SCAN_DIR, paste0(SCAN, ".", suffix, ".", chr, ".txt"))
    cat("Reading", f, "\n")
    read.table(f, header = TRUE)
  })
}

hap <- read_chr_files("scan")

snp_files <- file.path(SCAN_DIR, paste0(SCAN, ".snp_scan.", CHRS, ".txt"))
have_snp  <- all(file.exists(snp_files))
if (have_snp) {
  snp <- read_chr_files("snp_scan")
} else {
  cat("SNP scan files not yet available — skipping SNP panel\n")
  snp <- NULL
}

# ── Cumulative x-axis positions ───────────────────────────────────────────────
chr_order  <- CHRS
chr_colors <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
  chr_order
)
GAP_MB <- 5

chr_offsets <- bind_rows(
    hap %>% select(chr, pos),
    snp %>% select(chr, pos)
  ) %>%
  group_by(chr) %>%
  summarize(max_mb = max(pos) / 1e6, .groups = "drop") %>%
  filter(chr %in% chr_order) %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  arrange(chr) %>%
  mutate(offset = lag(cumsum(max_mb + GAP_MB), default = 0),
         mid    = offset + max_mb / 2)

add_xpos <- function(df) {
  df %>%
    filter(chr %in% chr_order) %>%
    mutate(chr = factor(chr, levels = chr_order)) %>%
    left_join(chr_offsets %>% select(chr, offset), by = "chr") %>%
    mutate(xpos = pos / 1e6 + offset)
}

hap <- add_xpos(hap)
if (have_snp) snp <- add_xpos(snp)

# ── Panel builder ─────────────────────────────────────────────────────────────
make_panel <- function(df, yvar, ylab, threshold = NULL) {
  p <- ggplot(df, aes(x = xpos, y = .data[[yvar]], color = chr)) +
    geom_point(size = 0.25, alpha = 0.5) +
    scale_color_manual(values = chr_colors, guide = "none") +
    scale_x_continuous(
      breaks = chr_offsets$mid,
      labels = chr_offsets$chr,
      expand = expansion(mult = 0.01)
    ) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(size = 9))
  if (!is.null(threshold))
    p <- p + geom_hline(yintercept = threshold, linetype = "dashed",
                        color = "grey50", linewidth = 0.4)
  p
}

# ── Build figure ──────────────────────────────────────────────────────────────
p1 <- make_panel(hap, "Wald_log10p", "Haplotype Wald  –log₁₀(p)", threshold = 10)
p3 <- make_panel(hap, "Falc_H2",     "Falconer H²")
p4 <- make_panel(hap, "Cutl_H2",     "Cutler H²")

if (have_snp) {
  p2  <- make_panel(snp, "Wald_log10p", "SNP Wald  –log₁₀(p)", threshold = 10)
  fig <- (p1 / p2 / p3 / p4)
} else {
  fig <- (p1 / p3 / p4)
}

fig <- fig +
  plot_annotation(
    title = "Malathion resistance — freqsmooth pipeline test (125 kb smoothing)",
    theme = theme(plot.title = element_text(size = 11, face = "bold"))
  )

# ── Save ──────────────────────────────────────────────────────────────────────
cat("Writing", OUT_FILE, "\n")
ggsave(OUT_FILE, fig, width = 8, height = 11, dpi = 150)
cat("Done.\n")
