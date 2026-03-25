SCAN_FILES <- c(
  "process/XQTL_talk/AGING/AGE_SY20_M_freqs250/AGE_SY20_M_freqs250.pseudoscan.txt"
)
SCAN_LABELS  <- c("Males day 20 selected vs control")
SCAN_COLOURS <- c("black")
THRESHOLD    <- 10
OUT_FILE     <- "figures/slides/aging.png"
FORMAT       <- "powerpoint"
OUT_HEIGHT_IN <- 5.8
PEAKS        <- NULL
GENES        <- NULL

source("plot_pseudoscan.R")
