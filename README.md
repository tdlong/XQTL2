# XQTL Pipeline

## Overview

This pipeline takes pooled-sequencing XQTL data from raw reads to genome scans and
summary figures. Two scan options are available:

- **Legacy scan** (`haps2scan.Apr2025`) — raw scan, no smoothing
- **Smooth scan** (`haps2scan.freqsmooth`) — stabilized scan with covariance and
  frequency smoothing (recommended for final analyses)

Both share the same upstream steps (alignment, REFALT counts, haplotype calling)
and the same downstream concat/plot step.

All scripts are run from the project root directory on a SLURM cluster.

---

## Step 1 — Get raw reads

The sequencing core will send download links. Save them to a file and generate a
download script:

```bash
cat links.txt | grep http | cut -f1 -d' ' | awk '{printf("wget %s\n",$0)}' > get_data.sh
```

Add a SLURM header to `get_data.sh` and submit. Store raw reads under `data/raw/<project_name>/`.

---

## Step 2 — Align reads to reference genome (fq → bam)

### Reference genome

Reference genome files go in `ref/` (too large for git). The pipeline expects `ref/dm6.fa`
with standard BWA and samtools indices. Copy from a shared location or index your own.

### Barcode-to-sample mapping file

Create a tab-delimited file mapping sequencing barcodes to sample names. Each row is one
sample with three fields: forward barcode, reverse barcode, sample name. Sample names
become the bam file prefixes and readgroup IDs used throughout the pipeline — choose them
carefully and consistently.

```
TGGCTATG    TTGTCAGC    R3con
GTCCTAGA    TTGTCAGC    R3age
ACTTGCCA    TTGTCAGC    R5con
TCTTCGTG    TTGTCAGC    R5age
```

Save this file to `helpfiles/` (e.g. `helpfiles/<project_name>/<project_name>.barcodes.txt`).

### Run alignment

```bash
mkdir -p data/bam/<project_name>
NN=$(wc -l < helpfiles/<project_name>/<project_name>.barcodes.txt)
sbatch --array=1-$NN scripts/fq2bam.sh \
    helpfiles/<project_name>/<project_name>.barcodes.txt \
    data/raw/<project_name> \
    data/bam/<project_name>
```

Bam files below ~1 GB likely indicate a failed library prep and should be reprocessed.

---

## Step 3 — Generate REFALT counts (bam → REFALT)

Create a file listing all bam paths for your experiment (pooled samples + founders).
Founders are pre-aligned; paths to the shared founder bams are in `helpfiles/founder.bams.txt`.

```bash
mkdir -p process/<project_name>
find data/bam/<project_name> -name "*.bam" -size +1G > helpfiles/<project_name>/bams
cat helpfiles/founder.bams.txt >> helpfiles/<project_name>/bams

sbatch scripts/bam2bcf2REFALT.sh \
    helpfiles/<project_name>/bams \
    process/<project_name>
```

---

## Step 4 — Call haplotypes (REFALT → haps)

### Haplotype parameters file

Create a parameter file for your experiment (e.g. `helpfiles/<project_name>/hap_params.R`).
This is an R script that is `source()`d by the pipeline. Required variables:

```r
# Founder set for this population
founders <- c("B1","B2","B3","B4","B5","B6","B7","AB8")

# Sample names — must exactly match the bam file prefixes / readgroup IDs from Step 2
names_in_bam <- c("R1con","R1age","R2con","R2age","R3con","R3age",
                   "R4con","R4age","R5con","R5age","R6con","R6age")

# Window step size in bp — haplotypes are inferred every step bp along the genome.
# 5000 (5 kb) is typical; use 10000 (10 kb) for very large experiments to save runtime.
step <- 5000

# Base half-window size in bp for haplotype inference.
# The caller uses an adaptive window: this value is the half-window at the
# chromosomal peak recombination rate. In low-recombination regions the window
# grows automatically (proportional to max_RR / local_RR) so that each window
# captures a similar number of informative recombination events regardless of
# local recombination rate. The polynomials describing recombination rate across
# dm6 chromosomes are hardcoded in REFALT2haps.code.R.
# 50000 (50 kb) is a reasonable base; windows in pericentromeric regions will
# typically be 5–10× larger.
size <- 50000

# Tree height cutoff: founders closer than this (Euclidean distance across SNPs in
# the window) are treated as indistinguishable. 2.5 is a conservative default.
h_cutoff <- 2.5
```

To generate `names_in_bam` from your bam directory:
```bash
echo -n "names_in_bam <- c(" && \
find data/bam/<project_name> -name "*.bam" -size +1G -print0 | \
xargs -0 -n1 basename | sed 's/.bam//' | sort | \
sed 's/.*/"&"/' | tr '\n' ',' | sed 's/,$//' && echo ")"
```

### Run haplotype calling

```bash
sbatch --array=1-5 scripts/REFALT2haps.sh \
    --parfile helpfiles/<project_name>/hap_params.R \
    --dir     process/<project_name>
```

---

## Step 5 — Run the scan

### Design file

Both scan options take a design file: a plain text table readable by R's `read.table()`,
saved with `write.table()` from R. Required columns:

| Column | Description |
|--------|-------------|
| `bam` | Sample name — must exactly match bam prefix / readgroup from Step 2 |
| `TRT` | `C` = control, `Z` = selected (rows with other values are ignored) |
| `REP` | Replicate number (integer) |
| `REPrep` | Technical replicate within replicate — usually `1` |
| `Num` | Number of flies in pool |
| `Proportion` | Proportion of flies selected (`NA` for controls) |

Extra columns are allowed and ignored. `Proportion` should be a decimal (not percent).
Create and save this file from R:

```r
design <- data.frame(
    bam        = c("R1con","R1age","R2con","R2age","R3con","R3age"),
    TRT        = c("C","Z","C","Z","C","Z"),
    REP        = c(1,1,2,2,3,3),
    REPrep     = 1,
    Num        = c(1205,115,1387,296,1631,174),
    Proportion = c(NA,0.087,NA,0.154,NA,0.088)
)
write.table(design, "helpfiles/<project_name>/design.txt")
```

### Legacy scan (no smoothing)

```bash
sbatch --array=1-5 scripts/haps2scan.Apr2025.sh \
    --rfile  helpfiles/<project_name>/design.txt \
    --dir    process/<project_name> \
    --outdir <scan_name>
```

### Smooth scan (recommended)

```bash
sbatch --array=1-5 scripts/haps2scan.freqsmooth.sh \
    --rfile          helpfiles/<project_name>/design.txt \
    --dir            process/<project_name> \
    --outdir         <scan_name> \
    --cov-smooth-kb  125 \
    --freq-smooth-kb 125
```

`--cov-smooth-kb` and `--freq-smooth-kb` specify the smoothing half-window in kilobases.
Set to `0` to disable that smoothing type. 125 kb is recommended for both.
The smooth scan also applies eigenvalue regularization to stabilize the Wald test.
Window counts are derived automatically from the actual data spacing in the haps file —
you always specify distances in kb.

The output directory `<scan_name>` is created inside `process/<project_name>/`. Choose a
name that reflects the analysis (e.g. `ZINC2_F_smooth125`).

### Per-chromosome outputs

Both pipelines write the same file layout — one pair of files per chromosome:

```
process/<project_name>/<scan_name>/
    <scan_name>.scan.chrX.txt
    <scan_name>.scan.chr2L.txt
    <scan_name>.scan.chr2R.txt
    <scan_name>.scan.chr3L.txt
    <scan_name>.scan.chr3R.txt
    <scan_name>.meansBySample.chrX.txt
    <scan_name>.meansBySample.chr2L.txt
    <scan_name>.meansBySample.chr2R.txt
    <scan_name>.meansBySample.chr3L.txt
    <scan_name>.meansBySample.chr3R.txt
```

If you run both pipelines on the same experiment (e.g. `ZINC2_F_legacy` and
`ZINC2_F_smooth125`) they each get their own subdirectory and are concatenated
independently in Step 6.

---

## Step 6 — Concatenate chromosomes and generate summary figures

Once all five chromosome jobs finish, pass the scan directory as the sole argument:

```bash
bash scripts/concat_Chromosome_Scans.sh process/<project_name>/<scan_name>
```

This merges the per-chromosome scan files, applies a light sliding-window smooth to
the summary statistics, and generates three figures:

| File | Description |
|------|-------------|
| `<scan_name>.scan.txt` | Full genome scan table |
| `<scan_name>.meansBySample.txt` | Per-founder frequency table |
| `<scan_name>.5panel.Mb.png` | 5-panel Manhattan plot (physical position) |
| `<scan_name>.5panel.cM.png` | 5-panel Manhattan plot (genetic position) |
| `<scan_name>.Manhattan.png` | Combined Manhattan plot |
| `<scan_name>.tar.gz` | All of the above bundled |

---

## Step 7 — Generate publication figures (smooth scan only)

The smooth scan has a dedicated plotting script (`scripts/plot_pseudoscan.R`) that
produces cleaner figures suitable for presentations and manuscripts. Each figure is
driven by a small R script stored in `helpfiles/<project_name>/`. Name the script
after the figure it produces — the `.R` and `.png` share a base name.

### Create a figure script

`helpfiles/<project_name>/<figure_name>.R` — single scan example:

```r
SCAN_FILES   <- c("process/<project_name>/<scan_name>/<scan_name>.scan.txt")
SCAN_COLOURS <- c("#1F78B4")
SCAN_LABELS  <- NULL       # NULL uses the scan file basename as the label
THRESHOLD    <- 10         # dashed line at this Wald -log10(p)
OUT_FILE     <- "process/<project_name>/<scan_name>/<figure_name>.png"
FORMAT       <- "powerpoint"
PEAKS        <- NULL
GENES        <- NULL

source("scripts/plot_pseudoscan.R")
```

`helpfiles/<project_name>/<figure_name>.R` — two-scan overlay (e.g. males and females):

```r
SCAN_FILES <- c(
    "process/<project_name>/<scan_name_M>/<scan_name_M>.scan.txt",
    "process/<project_name>/<scan_name_F>/<scan_name_F>.scan.txt"
)
SCAN_COLOURS <- c("#1F78B4", "#E31A1C")
SCAN_LABELS  <- c("Male", "Female")
THRESHOLD    <- 10
OUT_FILE     <- "process/<project_name>/<figure_name>.png"
FORMAT       <- "powerpoint"
PEAKS        <- NULL
GENES        <- NULL

source("scripts/plot_pseudoscan.R")
```

Height is calculated automatically (1.4 in per chromosome). Override with
`OUT_HEIGHT_IN <- 9.0` if needed.

Gene and peak annotations both use Mb coordinates:

```r
GENES <- data.frame(
    name   = c("Ace",   "Cyp6g1"),
    chr    = c("chr3R", "chr2R"),
    pos_mb = c(9.07,    12.19)
)
PEAKS <- data.frame(
    label  = c("peak1"),
    chr    = c("chr3R"),
    pos_mb = c(9.1)
)
```

### FORMAT options

| FORMAT | Width | DPI | Use for |
|--------|-------|-----|---------|
| `manuscript_half` | 3.5 in | 300 | half-width journal figure |
| `manuscript_full` | 7.0 in | 300 | full-width journal figure |
| `manuscript_half_hires` | 3.5 in | 600 | high-res submission |
| `manuscript_full_hires` | 7.0 in | 600 | high-res submission |
| `powerpoint` | 8.0 in | 150 | slides |
| `web` | 7.0 in | 150 | web/HTML |
| `email` | 6.0 in | 100 | email preview |

### Run the figure script

```bash
Rscript helpfiles/<project_name>/<figure_name>.R
```

---

## Step 8 — Download and explore results

The concat step (Step 6) bundles the scan table, means table, and summary figures
into a single tarball. Download it in one shot:

```bash
scp <user>@<cluster>:<project_path>/process/<project_name>/<scan_name>/<scan_name>.tar.gz .
tar -xzf <scan_name>.tar.gz
```

The tarball contains:

| File | Description |
|------|-------------|
| `<scan_name>.scan.txt` | Full genome scan table |
| `<scan_name>.meansBySample.txt` | Per-founder frequency table |
| `<scan_name>.5panel.Mb.png` | 5-panel Manhattan (physical position) |
| `<scan_name>.5panel.cM.png` | 5-panel Manhattan (genetic position) |
| `<scan_name>.Manhattan.png` | Combined Manhattan plot |

---

## Interactive plotting functions

`scripts/XQTL_plotting_functions.R` provides functions for interactive exploration
of scan results (works with both pipeline outputs). Run from the project root or
point the paths at local copies of the files downloaded in Step 8.

```r
library(tidyverse)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

source("scripts/XQTL_plotting_functions.R")

scan_dir <- "process/<project_name>/<scan_name>"
df1 <- as_tibble(read.table(file.path(scan_dir, "<scan_name>.scan.txt")))
df2 <- as_tibble(read.table(file.path(scan_dir, "<scan_name>.meansBySample.txt")))

# Genome-wide Manhattan plots
XQTL_Manhattan_5panel(df1, cM = FALSE)
XQTL_Manhattan_5panel(df1, cM = TRUE)
XQTL_Manhattan(df1, cM = FALSE, color_scheme = "UCI")

# Regional plots
XQTL_region(df1, "chr3R", 18250000, 19000000, "Wald_log10p")
XQTL_change_average(df2, "chr3R", 18250000, 19000000)
XQTL_change_average(df2, "chr3R", 18250000, 19000000, reference_strain = "B5")
XQTL_change_byRep(df2, "chr3R", 18250000, 19000000)
XQTL_beforeAfter_selectReps(df2, "chr3R", 18250000, 19000000, reps = c(1,7,9,12))
XQTL_combined_plot(df1, df2, "chr3R", 18250000, 19000000)

# Zoom to a peak automatically — adjust drop thresholds until the window looks right
out <- XQTL_zoom(df1, "chr2L", 15000000, 16000000, drop_left = 3, drop_right = 3)
out$plot
A1 <- XQTL_region(df1, out$chr, out$start, out$stop, "Wald_log10p")
A2 <- XQTL_change_average(df2, out$chr, out$start, out$stop)
A1 / A2

# With gene track (requires rtracklayer / GenomicRanges and a GTF file)
# wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/genes/dm6.ncbiRefSeq.gtf.gz
gtf <- rtracklayer::import("dm6.ncbiRefSeq.gtf")
A3  <- XQTL_genes(gtf, out$chr, out$start, out$stop)
A1 / A2 / A3
```

---

## Putting it all together

For a new experiment, copy the block below to
`scripts_oneoffs/<project_name>_pipeline.sh`, fill in the variables at the top,
and run it from the project root. Each step is submitted as a SLURM job with
`afterok` dependencies so the whole pipeline runs unattended.

```bash
#!/bin/bash
# Edit these variables, then: bash scripts_oneoffs/<project_name>_pipeline.sh

PROJECT=<project_name>
BARCODES=helpfiles/${PROJECT}/${PROJECT}.barcodes.txt
PARFILE=helpfiles/${PROJECT}/hap_params.R
COV_KB=125
FREQ_KB=125

# One entry per scan (sex, treatment, etc.) — add or remove as needed
DESIGNS=(  helpfiles/${PROJECT}/design_male.txt   helpfiles/${PROJECT}/design_female.txt )
OUTDIRS=(  ${PROJECT}_M_smooth${COV_KB}           ${PROJECT}_F_smooth${COV_KB}           )

# Figure scripts that run after all scans are concatenated
FIGURES=(  helpfiles/${PROJECT}/${PROJECT}_MF_smooth${COV_KB}.R )

# ── Step 2: align reads ───────────────────────────────────────────────────────
NN=$(wc -l < ${BARCODES})
mkdir -p data/bam/${PROJECT}
jid_bam=$(sbatch --parsable --array=1-${NN} scripts/fq2bam.sh \
    ${BARCODES} data/raw/${PROJECT} data/bam/${PROJECT})
echo "fq2bam:      $jid_bam"

# ── Step 3: REFALT counts ─────────────────────────────────────────────────────
mkdir -p process/${PROJECT}
find data/bam/${PROJECT} -name "*.bam" -size +1G > helpfiles/${PROJECT}/bams
cat helpfiles/founder.bams.txt >> helpfiles/${PROJECT}/bams

jid_refalt=$(sbatch --parsable --dependency=afterok:${jid_bam} \
    --array=1-5 scripts/bam2bcf2REFALT.sh \
    helpfiles/${PROJECT}/bams process/${PROJECT})
echo "REFALT:      $jid_refalt"

# ── Step 4: haplotypes ────────────────────────────────────────────────────────
jid_haps=$(sbatch --parsable --dependency=afterok:${jid_refalt} \
    --array=1-5 scripts/REFALT2haps.sh \
    --parfile ${PARFILE} --dir process/${PROJECT})
echo "haps:        $jid_haps"

# ── Steps 5-6: scans + concat (one pair per design, all share same haps) ─────
jid_concats=""
for i in "${!DESIGNS[@]}"; do
    design=${DESIGNS[$i]}
    outdir=${OUTDIRS[$i]}

    jid_scan=$(sbatch --parsable --dependency=afterok:${jid_haps} \
        --array=1-5 scripts/haps2scan.freqsmooth.sh \
        --rfile ${design} --dir process/${PROJECT} --outdir ${outdir} \
        --cov-smooth-kb ${COV_KB} --freq-smooth-kb ${FREQ_KB})
    echo "scan ${outdir}: $jid_scan"

    jid_concat=$(sbatch --parsable --dependency=afterok:${jid_scan} \
        -A tdlong_lab -p standard --mem=10G \
        --wrap="bash scripts/concat_Chromosome_Scans.sh process/${PROJECT}/${outdir}")
    echo "concat ${outdir}: $jid_concat"

    jid_concats="${jid_concats}:${jid_concat}"
done

# ── Step 7: figures (after all concats finish) ────────────────────────────────
dep="afterok${jid_concats}"
for fig in "${FIGURES[@]}"; do
    jid_fig=$(sbatch --parsable --dependency=${dep} \
        -A tdlong_lab -p standard --mem=10G \
        --wrap="module load R/4.2.2 && Rscript ${fig}")
    echo "figure ${fig}: $jid_fig"
done
```

---

## Directory structure

```
XQTL2/
├── scripts/              # Core pipeline scripts (tracked in git)
├── scripts_oneoffs/      # Experiment-specific submit scripts (not tracked)
├── helpfiles/
│   ├── flymap.r6.txt                         (tracked)
│   ├── founder.bams.txt                      (tracked)
│   └── <project_name>/
│       ├── <project_name>.barcodes.txt       (Step 2)
│       ├── bams                              (Step 3)
│       ├── hap_params.R                      (Step 4)
│       ├── design.txt                        (Step 5)
│       └── <figure_name>.R                   (Step 7, one per figure)
├── data/
│   ├── raw/<project_name>/                   (Step 1 — raw reads)
│   └── bam/<project_name>/                   (Step 2 — aligned bams)
├── ref/                  # Reference genome (not tracked)
├── process/
│   └── <project_name>/
│       ├── RefAlt.<chr>.txt                  (Step 3 — allele counts)
│       ├── calls.<chr>.bcf                   (Step 3 — intermediate)
│       ├── R.haps.<chr>.rds                  (Step 4 — SNP table)
│       ├── R.haps.<chr>.out.rds              (Step 4 — haplotype estimates)
│       └── <scan_name>/                      (one per scan run, Steps 5-6)
│           ├── <scan_name>.scan.<chr>.txt
│           ├── <scan_name>.meansBySample.<chr>.txt
│           ├── <scan_name>.scan.txt          (after concat)
│           ├── <scan_name>.meansBySample.txt (after concat)
│           ├── <scan_name>.5panel.Mb.png
│           ├── <scan_name>.5panel.cM.png
│           ├── <scan_name>.Manhattan.png
│           └── <scan_name>.tar.gz
└── figures/              # Publication figures from configs/*.R (not tracked)
```

To check what exists for a given project on the cluster:

```bash
bash scripts/show_project_layout.sh <project_name>
```
