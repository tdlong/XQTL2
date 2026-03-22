###############################################################################
# plot_H2_overlay.R  —  source this file, do not edit it
#
# Overlays Falconer H² and Cutler H² from the same scan file on one figure.
# Same 5-panel per-chromosome layout as plot_freqsmooth_H2.R.
#
# Required:
#   SCAN_FILE     single .scan.txt or .snp_scan.txt path
#   OUT_FILE      output .png path
#   FORMAT        powerpoint | manuscript_half | manuscript_full | web | email
#
# Optional:
#   THRESHOLD     numeric — dashed line (or NULL)
#   OUT_HEIGHT_IN numeric — overrides 1.4 in per chromosome default
#   HET_BOUNDS    data frame with columns chr, eu_start, eu_end
#   FALC_COLOUR   colour for Falconer H² (default "#E31A1C")
#   CUTL_COLOUR   colour for Cutler H²   (default "#33A02C")
###############################################################################

library(tidyverse)

FORMAT_PRESETS <- list(
  manuscript_half       = c(w = 3.5, dpi = 300, font =  7),
  manuscript_full       = c(w = 7.0, dpi = 300, font =  8),
  manuscript_half_hires = c(w = 3.5, dpi = 600, font =  7),
  manuscript_full_hires = c(w = 7.0, dpi = 600, font =  8),
  powerpoint            = c(w = 8.0, dpi = 150, font =  9),
  web                   = c(w = 7.0, dpi = 150, font = 10),
  email                 = c(w = 6.0, dpi = 100, font =  9)
)

if (!FORMAT %in% names(FORMAT_PRESETS))
  stop("Unknown FORMAT '", FORMAT, "'. Choose from: ",
       paste(names(FORMAT_PRESETS), collapse = ", "))

fmt       <- FORMAT_PRESETS[[FORMAT]]
BASE_FONT <- fmt["font"]

chr_order  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
chr_labels <- c(chrX = "X", chr2L = "2L", chr2R = "2R", chr3L = "3L", chr3R = "3R")

if (!exists("THRESHOLD"))    THRESHOLD    <- NULL
if (!exists("FALC_COLOUR"))  FALC_COLOUR  <- "#E31A1C"
if (!exists("CUTL_COLOUR"))  CUTL_COLOUR  <- "#33A02C"

if (!exists("HET_BOUNDS")) {
  HET_BOUNDS <- tribble(
    ~chr,    ~eu_start, ~eu_end,
    "chrX",   2.5,      21.2,
    "chr2L",  0.5,      22.9,
    "chr2R",  1.3,      25.1,
    "chr3L",  0.7,      24.0,
    "chr3R",  4.5,      32.0
  )
}

if (!file.exists(SCAN_FILE)) stop("File not found: ", SCAN_FILE)

scan_df <- read.table(SCAN_FILE, header = TRUE) %>%
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

if (nrow(scan_df) == 0) stop("No H2 data loaded — check SCAN_FILE path.")

if (!exists("OUT_HEIGHT_IN"))
  OUT_HEIGHT_IN <- length(unique(as.character(scan_df$chr))) * 1.4

colour_map <- c("Falconer H\u00b2" = FALC_COLOUR, "Cutler H\u00b2" = CUTL_COLOUR)

chr_label_df <- scan_df %>%
  distinct(chr) %>%
  mutate(label_text = chr_labels[as.character(chr)])

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

if (!is.null(THRESHOLD))
  p <- p + geom_hline(yintercept = THRESHOLD, linetype = "dashed",
                      colour = "grey40", linewidth = 0.4)

png(OUT_FILE, width = fmt["w"], height = OUT_HEIGHT_IN,
    units = "in", res = fmt["dpi"])
print(p)
dev.off()
cat("Saved:", OUT_FILE, "\n")
