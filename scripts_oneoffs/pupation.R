SCAN_FILES <- c(
  "process/XQTL_talk/PUPATION/PUPATION_TB_freqs250/PUPATION_TB_freqs250.pseudoscan.txt"
)
SCAN_LABELS  <- c("Top vs Bottom")
SCAN_COLOURS <- c("black")
THRESHOLD    <- 10
OUT_FILE     <- "figures/slides/pupation.png"
FORMAT       <- "powerpoint"
OUT_HEIGHT_IN <- 5.8
PEAKS        <- NULL
GENES        <- NULL

source("plot_pseudoscan.R")
