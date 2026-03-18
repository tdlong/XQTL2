# XQTL Pipeline

## Overview

This pipeline takes pooled-sequencing XQTL data from raw reads to genome scans and summary figures. Two scan pipelines are available:

- **Legacy scan** (`haps2scan.Apr2025`) — raw scan, no smoothing
- **Smooth scan** (`haps2scan.freqsmooth`) — stabilized scan with covariance and frequency smoothing (recommended for final analyses)

Both share the same upstream steps (alignment, REFALT counts, haplotype calling) and the same downstream concat/plot step.

---

## Step 1 — Get raw reads

```bash
# Save the email from the sequencing core as blah.txt and extract download links
cat blah.txt | grep http | cut -f1 -d' ' | awk '{printf("wget %s\n",$0)}' > get_data.sh
mkdir data/raw/Oct28_24
# Add SLURM header to get_data.sh, then:
sbatch get_data.sh
```

---

## Step 2 — Align reads to reference genome (fq → bam)

```bash
# Reference genome files go in ref/ (too large to track in git)
# Copy from shared location or download dm6
cp /dfs7/adl/tdlong/fly_pool/newpipeline_aging/XQTL_pipeline/ref/* ref/.

# Create a barcode-to-sample mapping file (tab-delimited: F_barcode R_barcode sample_name)
# e.g. helperfiles/readname.mapping.Oct28.txt

mkdir data/bam/Oct28_24
NN=`wc -l helperfiles/readname.mapping.Oct28.txt | cut -f1 -d' '`
sbatch --array=1-$NN scripts/fq2bam.sh \
    helperfiles/readname.mapping.Oct28.txt data/raw/Oct28_24 data/bam/Oct28_24
```

Bam files below ~1G likely indicate failed library prep and should be reprocessed.

---

## Step 3 — Generate REFALT counts (bam → REFALT)

```bash
mkdir process/Oct28_24
find data/bam/Oct28_24 -name "*.bam" -size +1G > helpfiles/Oct28_24.bams
cat helpfiles/founder.bams.txt | grep "B" >> helpfiles/Oct28_24.bams
sbatch scripts/bam2bcf2REFALT.sh helpfiles/Oct28_24.bams process/Oct28_24
```

---

## Step 4 — Call haplotypes (REFALT → haps)

Edit `helpfiles/haplotype_parameters.R` to reflect your experiment (founders, sample names, step size, window size, tree cutoff). Then:

```bash
sbatch scripts/REFALT2haps.Andreas.sh helpfiles/haplotype_parameters.R process/Oct28_24
```

Key parameters in `haplotype_parameters.R`:

```r
founders  <- c("B1","B2","B3","B4","B5","B6","B7","AB8")
names_in_bam <- c("R1con","R1age", ...)   # must match bam file prefixes
step      <- 5000    # window step in bp (5 kb typical; 10 kb for large experiments)
size      <- 50000   # half-window size in bp for haplotype inference
h_cutoff  <- 2.5     # tree height cutoff to flag indistinguishable founders
```

---

## Step 5 — Run the scan

### Design file

Both scan pipelines take a design file readable by `read.table()` with these required columns:

| Column | Description |
|--------|-------------|
| `bam` | Sample name (must match bam file prefix / readgroup) |
| `TRT` | `C` = control, `Z` = selected (other values ignored) |
| `REP` | Replicate number |
| `REPrep` | Technical replicate within replicate (usually `1`) |
| `Num` | Number of flies in pool |
| `Proportion` | Proportion selected (`NA` for controls) |

Example:
```
bam         TRT  REP  REPrep  Num   Proportion
STV1_F_Con  C    1    1       1205  NA
STV1_F_Res  Z    1    1       115   0.0871
STV2_F_Con  C    2    1       1387  NA
STV2_F_Res  Z    2    1       296   0.1540
```

### Legacy scan (no smoothing)

```bash
sbatch --array=1-5 scripts/haps2scan.Apr2025.sh \
    helpfiles/mydesign.txt process/Oct28_24 SCAN_NAME
```

### Smooth scan (recommended)

```bash
# COV_SMOOTH_KB:  covariance smoothing half-window in kb (0 = off; 125 recommended)
# FREQ_SMOOTH_KB: frequency smoothing half-window in kb  (0 = off; 125 recommended)
sbatch --array=1-5 scripts/haps2scan.freqsmooth.sh \
    helpfiles/mydesign.txt process/Oct28_24 SCAN_NAME 125 125
```

Pipeline 2 applies:
1. **Covariance smoothing** — running mean of reconstruction error covariance matrices
   across a ±`COV_SMOOTH_KB` half-window. Reduces noise in the Wald test statistic.
2. **Eigenvalue regularization** — floors small eigenvalues to stabilize matrix inversion.
3. **Frequency smoothing** (optional) — running mean of haplotype frequency vectors
   across a ±`FREQ_SMOOTH_KB` half-window before computing the Wald test.

Window counts are derived automatically from the actual data spacing — you always
specify distances in kb.

See `scripts_oneoffs/run_all_scans_125.sh` for an example of how to chain multiple
experiments with dependency management.

---

## Step 6 — Concatenate chromosomes and generate summary figures

```bash
bash scripts/concat_Chromosome_Scans.sh process/Oct28_24/SCAN_NAME
```

This merges the per-chromosome scan files into a single table, applies a light
sliding-window smooth to the summary statistics, and generates three figures:

- `SCAN_NAME.5panel.Mb.png` — 5-panel Manhattan plot (physical position)
- `SCAN_NAME.5panel.cM.png` — 5-panel Manhattan plot (genetic position)
- `SCAN_NAME.Manhattan.png` — combined Manhattan plot

Output files are also bundled into `SCAN_NAME.tar.gz`.

---

## Step 7 — Download results and explore

```bash
scp tdlong@hpc3.rcic.uci.edu:.../SCAN_NAME/SCAN_NAME.scan.txt .
scp tdlong@hpc3.rcic.uci.edu:.../SCAN_NAME/SCAN_NAME.meansBySample.txt .
scp tdlong@hpc3.rcic.uci.edu:.../SCAN_NAME/SCAN_NAME.5panel.cM.png .
```

The `.scan.txt` and `.meansBySample.txt` files are the primary outputs for downstream analysis.

---

## Plotting functions

```r
library(tidyverse)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

source("scripts/XQTL_plotting_functions.R")
df1 <- as_tibble(read.table("SCAN_NAME.scan.txt"))
df2 <- as_tibble(read.table("SCAN_NAME.meansBySample.txt"))

XQTL_Manhattan_5panel(df1, cM = FALSE)
XQTL_Manhattan_5panel(df1, cM = TRUE)
XQTL_Manhattan(df1, cM = FALSE, color_scheme = "UCI")
XQTL_change_average(df2, "chr3R", 18250000, 19000000)
XQTL_change_average(df2, "chr3R", 18250000, 19000000, reference_strain = "B5")
XQTL_change_byRep(df2, "chr3R", 18250000, 19000000)
XQTL_beforeAfter_selectReps(df2, "chr3R", 18250000, 19000000, reps = c(1,7,9,12))
XQTL_region(df1, "chr3R", 18250000, 19000000, "Wald_log10p")
XQTL_combined_plot(df1, df2, "chr3R", 18250000, 19000000)

# Zoom to a peak automatically
out <- XQTL_zoom(df1, "chr2L", 15000000, 16000000, 3, 3)
out$plot
```

---

## Directory structure

```
XQTL2/
├── scripts/          # Core pipeline scripts (tracked)
├── scripts_oneoffs/  # Experiment-specific submit and plot scripts (not tracked)
├── helpfiles/        # Design files, bam lists, parameter files
├── configs/          # Experiment-specific plotting configs
├── data/             # Raw and aligned data (not tracked)
├── ref/              # Reference genome (not tracked)
├── process/          # Pipeline outputs (not tracked)
└── figures/          # Summary figures (not tracked)
```
