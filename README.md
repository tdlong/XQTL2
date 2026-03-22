# XQTL Pipeline

## Overview

This pipeline takes pooled-sequencing XQTL data from raw reads to genome scans and
publication figures. All scripts are run from the project root on a SLURM cluster.

**Pipeline at a glance:**

1. Get raw reads
2. Align reads (`fq2bam.sh`)
3. Generate allele counts (`bam2bcf2REFALT.sh`)
4. Call haplotypes (`REFALT2haps.sh`)
5. Run the scan
   - 5a. Haplotype scan (`run_scan.sh` — one command does smoothing, scanning, and concatenation)
   - 5b. SNP scan (`run_snp_scan.sh` — optional add-on after 5a)
6. Generate figures

A legacy scan without smoothing (`haps2scan.Apr2025.sh`) is also available.

---

## Step 1 — Get raw reads

The sequencing core will send download links. Save them to a file and generate a
download script:

```bash
cat links.txt | grep http | cut -f1 -d' ' | awk '{printf("wget %s\n",$0)}' > get_data.sh
```

Add a SLURM header to `get_data.sh` and submit. Store raw reads under `data/raw/<project>/`.

---

## Step 2 — Align reads (fq to bam)

### Reference genome

Reference genome files go in `ref/` (too large for git). The pipeline expects `ref/dm6.fa`
with standard BWA and samtools indices. Copy from a shared location or index your own.

### Barcode-to-sample mapping file

Create a tab-delimited file mapping sequencing barcodes to sample names. Each row is one
sample with three fields: forward barcode, reverse barcode, sample name. Sample names
become the bam file prefixes and readgroup IDs used throughout the pipeline.

```
TGGCTATG    TTGTCAGC    R3con
GTCCTAGA    TTGTCAGC    R3age
ACTTGCCA    TTGTCAGC    R5con
TCTTCGTG    TTGTCAGC    R5age
```

Save this file to `helpfiles/<project>/<project>.barcodes.txt`.

### Run alignment

```bash
mkdir -p data/bam/<project>
NN=$(wc -l < helpfiles/<project>/<project>.barcodes.txt)
sbatch --array=1-$NN scripts/fq2bam.sh \
    helpfiles/<project>/<project>.barcodes.txt \
    data/raw/<project> \
    data/bam/<project>
```

Bam files below ~1 GB likely indicate a failed library prep and should be reprocessed.

---

## Step 3 — Generate REFALT counts (bam to REFALT)

Create a file listing all bam paths for your experiment (pooled samples + founders).
Founders are pre-aligned; paths to the shared founder bams are in `helpfiles/founder.bams.txt`.

```bash
mkdir -p process/<project>
find data/bam/<project> -name "*.bam" -size +1G > helpfiles/<project>/bams
cat helpfiles/founder.bams.txt >> helpfiles/<project>/bams

sbatch scripts/bam2bcf2REFALT.sh \
    helpfiles/<project>/bams \
    process/<project>
```

---

## Step 4 — Call haplotypes (REFALT to haps)

### Haplotype parameters file

Create `helpfiles/<project>/hap_params.R`:

```r
# Founder set for this population
founders <- c("B1","B2","B3","B4","B5","B6","B7","AB8")

# Sample names — must exactly match bam prefixes from Step 2
names_in_bam <- c("R1con","R1age","R2con","R2age","R3con","R3age",
                   "R4con","R4age","R5con","R5age","R6con","R6age")

# Window step size in bp (5000 typical; 10000 for very large experiments)
step <- 5000

# Base half-window in bp for haplotype inference.
# The caller adapts this: in low-recombination regions the window grows
# proportional to max_RR / local_RR, so each window captures similar
# recombination events regardless of position.
size <- 50000

# Tree height cutoff for founder distinguishability (2.5 is default)
h_cutoff <- 2.5
```

To generate `names_in_bam` from your bam directory:
```bash
echo -n "names_in_bam <- c(" && \
find data/bam/<project> -name "*.bam" -size +1G -print0 | \
xargs -0 -n1 basename | sed 's/.bam//' | sort | \
sed 's/.*/"&"/' | tr '\n' ',' | sed 's/,$//' && echo ")"
```

### Run haplotype calling

```bash
sbatch --array=1-5 scripts/REFALT2haps.sh \
    --parfile helpfiles/<project>/hap_params.R \
    --dir     process/<project>
```

---

## Step 5a — Haplotype scan

`run_scan.sh` handles everything: smoothing haplotype frequencies, running the
Wald test and heritability estimates, concatenating chromosomes, and optionally
generating figures. One command per scan.

### Design file

Create a plain text table with one row per sample. Required columns:

| Column | Description |
|--------|-------------|
| `bam` | Sample name (must match bam prefix from Step 2) |
| `TRT` | `C` = control, `Z` = selected |
| `REP` | Replicate number (integer) |
| `REPrep` | Technical replicate within replicate (usually `1`) |
| `Num` | Number of flies in pool |
| `Proportion` | Fraction selected (`NA` for controls) |

Create and save from R:

```r
design <- data.frame(
    bam        = c("R1con","R1age","R2con","R2age","R3con","R3age"),
    TRT        = c("C","Z","C","Z","C","Z"),
    REP        = c(1,1,2,2,3,3),
    REPrep     = 1,
    Num        = c(1205,115,1387,296,1631,174),
    Proportion = c(NA,0.087,NA,0.154,NA,0.088)
)
write.table(design, "helpfiles/<project>/design.txt")
```

### Run the scan

```bash
bash scripts/run_scan.sh \
    --design    helpfiles/<project>/design.txt \
    --dir       process/<project> \
    --scan      <scan_name> \
    --after     $JID_HAPS
```

That's it. This submits all SLURM jobs (smooth, hap scan, concat) with proper
dependency chaining.

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--design` | (required) | Path to design file |
| `--dir` | (required) | Project directory (e.g. `process/<project>`) |
| `--scan` | (required) | Scan name — becomes output subdirectory |
| `--smooth` | 250 | Smoothing half-window in kb |
| `--figure` | (none) | R figure script to run after concat |
| `--after` | (none) | SLURM job ID to wait on before starting |

### Smoothing window

The default smoothing window is 250 kb. This was chosen by comparing 125 kb and
250 kb on real data (malathion experiment): plotting smoothed founder haplotype
frequencies versus genomic position and comparing against expectations from
simulations. At 250 kb the founder frequency estimates are stable without
over-smoothing genuine biological signal. You can override with `--smooth 125`
or any other value.

### What run_scan.sh submits

For reference, `run_scan.sh` chains these SLURM jobs automatically:

1. `smooth_haps.sh` — smooth haplotype frequencies and covariances (5 array tasks)
2. `hap_scan.sh` — Wald test + heritability at each haplotype window (5 array tasks, after #1)
3. `concat_scans.sh` — merge per-chromosome files and generate Manhattan plots (after #2)
4. Figure script (after #3; only if `--figure`)

### Legacy scan (alternative — no smoothing)

```bash
sbatch --array=1-5 scripts/haps2scan.Apr2025.sh \
    --rfile  helpfiles/<project>/design.txt \
    --dir    process/<project> \
    --outdir <scan_name>
```

Followed by `bash scripts/concat_scans.sh process/<project>/<scan_name>` to
concatenate chromosomes.

---

## Step 5b — SNP scan (optional)

The SNP scan imputes per-SNP allele frequencies from the smoothed haplotype
estimates produced in Step 5a, then runs a Wald test (df=1) at every SNP. This
gives much finer genomic resolution than the haplotype scan, at the cost of
relying on imputation from the founder haplotypes.

The SNP scan uses the same smoothed data as the haplotype scan, so `run_scan.sh`
(Step 5a) must have already run. Use `--after` to chain it after a running scan,
or run it any time after the haplotype scan has completed.

### SNP table

The scan requires a table of per-founder allele frequencies at every SNP. This is
a one-time preparation step per population — see `helpfiles/snp_tables/` for
details.

### Run the SNP scan

```bash
bash scripts/run_snp_scan.sh \
    --design    helpfiles/<project>/design.txt \
    --dir       process/<project> \
    --scan      <scan_name> \
    --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founders  A1,A2,A3,A4,A5,A6,A7,AB8 \
    --figure    helpfiles/<project>/snp_figure.R
```

Use the same `--scan` name as Step 5a — the SNP scan output goes into the same
directory alongside the haplotype scan results.

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--design` | (required) | Path to design file |
| `--dir` | (required) | Project directory |
| `--scan` | (required) | Scan name (same as Step 5a) |
| `--snp-table` | (required) | SNP frequency table |
| `--founders` | (required) | Comma-separated founder names matching SNP table columns |
| `--figure` | (none) | R figure script to run after concat |
| `--after` | (none) | SLURM job ID to wait on (e.g. from run_scan.sh) |

---

## Step 6 — Generate publication figures

Each figure is controlled by a short R parameter file that you write. The
parameter file sets variables (input paths, colors, output path, etc.) and then
calls `source()` on one of the plot engines in `scripts/`. You run the parameter
file with `Rscript` and it produces the figure.

### How it works

1. Write a parameter file, e.g. `helpfiles/<project>/wald_figure.R`
2. Run it: `Rscript helpfiles/<project>/wald_figure.R`
3. The plot engine reads the scan data, makes the figure, and saves the PNG

You can pass the parameter file to `run_scan.sh --figure` or
`run_snp_scan.sh --figure` so that the figure is generated automatically
as part of the pipeline.

### Plot engines

Three plot engines are available:

**`scripts/plot_pseudoscan.R`** — 5-panel Wald -log10(p) line plot from the haplotype scan.

**`scripts/plot_H2_overlay.R`** — 5-panel Falconer + Cutler heritability overlay from the haplotype scan.

**`scripts/plot_freqsmooth_snp.R`** — 5-panel Wald -log10(p) dot plot from the SNP scan.

### Example parameter files

**Haplotype Wald scan** — save as `helpfiles/<project>/wald_figure.R`:

```r
# Inputs
SCAN_FILES   <- c("process/<project>/<scan_name>/<scan_name>.scan.txt")
SCAN_COLOURS <- c("#1F78B4")
SCAN_LABELS  <- NULL       # NULL uses the file basename
THRESHOLD    <- 10         # dashed horizontal line at this -log10(p)

# Output
OUT_FILE     <- "process/<project>/<scan_name>/wald.png"
FORMAT       <- "powerpoint"    # sets width and DPI (see table below)

# Optional annotations (set to NULL to skip)
PEAKS        <- NULL
GENES        <- NULL

# Run the plot engine
source("scripts/plot_pseudoscan.R")
```

**Two-scan overlay** (e.g. male vs female) — save as `helpfiles/<project>/MF_overlay.R`:

```r
SCAN_FILES   <- c("process/<project>/<scan_M>/<scan_M>.scan.txt",
                   "process/<project>/<scan_F>/<scan_F>.scan.txt")
SCAN_COLOURS <- c("#1F78B4", "#E31A1C")
SCAN_LABELS  <- c("Male", "Female")
THRESHOLD    <- 10
OUT_FILE     <- "process/<project>/MF_overlay.png"
FORMAT       <- "powerpoint"
PEAKS        <- NULL
GENES        <- NULL

source("scripts/plot_pseudoscan.R")
```

**Heritability overlay** — save as `helpfiles/<project>/H2_figure.R`:

```r
SCAN_FILE <- "process/<project>/<scan_name>/<scan_name>.scan.txt"
OUT_FILE  <- "process/<project>/<scan_name>/H2.png"
FORMAT    <- "powerpoint"

source("scripts/plot_H2_overlay.R")
```

**SNP scan** — save as `helpfiles/<project>/snp_figure.R`:

```r
SCAN_FILE  <- "process/<project>/<scan_name>/<scan_name>.snp_scan.txt"
OUT_FILE   <- "process/<project>/<scan_name>/snp_wald.png"
FORMAT     <- "powerpoint"
THRESHOLD  <- 10

source("scripts/plot_freqsmooth_snp.R")
```

### Run a figure script

On the cluster:

```bash
module load R/4.2.2
Rscript helpfiles/<project>/wald_figure.R
```

Or pass it to the scan driver so it runs automatically after concat:

```bash
bash scripts/run_scan.sh \
    --design ... --dir ... --scan ... \
    --figure helpfiles/<project>/wald_figure.R
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

### Optional annotations

Gene and peak labels use Mb coordinates. Add these to any parameter file
that uses `plot_pseudoscan.R`:

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

---

## Step 7 — Download and explore results

The concat step bundles scan tables and quick-look figures into a tarball:

```bash
scp <user>@<cluster>:<project_path>/process/<project>/<scan_name>/<scan_name>.tar.gz .
tar -xzf <scan_name>.tar.gz
```

For interactive exploration, `scripts/XQTL_plotting_functions.R` provides functions
that work with both pipeline outputs:

```r
library(tidyverse)
library(patchwork)
source("scripts/XQTL_plotting_functions.R")

scan_dir <- "process/<project>/<scan_name>"
df1 <- as_tibble(read.table(file.path(scan_dir, "<scan_name>.scan.txt")))
df2 <- as_tibble(read.table(file.path(scan_dir, "<scan_name>.meansBySample.txt")))

# Genome-wide Manhattan
XQTL_Manhattan_5panel(df1, cM = FALSE)

# Regional plots
XQTL_region(df1, "chr3R", 18250000, 19000000, "Wald_log10p")
XQTL_change_average(df2, "chr3R", 18250000, 19000000)
XQTL_combined_plot(df1, df2, "chr3R", 18250000, 19000000)

# Zoom to a peak automatically
out <- XQTL_zoom(df1, "chr2L", 15000000, 16000000, drop_left = 3, drop_right = 3)
A1 <- XQTL_region(df1, out$chr, out$start, out$stop, "Wald_log10p")
A2 <- XQTL_change_average(df2, out$chr, out$start, out$stop)
A1 / A2
```

---

## Worked example — start-to-finish for one experiment

This example shows every step for a typical experiment with two conditions
(e.g. male and female). You submit the script once and it runs unattended on the
cluster. At the end you have, for each condition:

- A genome-wide haplotype scan with Wald -log10(p) and heritability estimates
- Manhattan plots
- Optionally, a per-SNP scan at single-nucleotide resolution
- Publication figures

Before running, you need to have prepared the input files described in
Steps 1–4 above: barcodes, haplotype parameters, design files, and figure
parameter files (Step 6).

Copy this to `scripts_oneoffs/<project>_pipeline.sh`, fill in the variables,
and run with `bash scripts_oneoffs/<project>_pipeline.sh`.

```bash
#!/bin/bash
PROJECT=<project>
BARCODES=helpfiles/${PROJECT}/${PROJECT}.barcodes.txt
PARFILE=helpfiles/${PROJECT}/hap_params.R

# One entry per condition
DESIGNS=(  helpfiles/${PROJECT}/design_male.txt   helpfiles/${PROJECT}/design_female.txt )
OUTDIRS=(  ${PROJECT}_M_smooth250                 ${PROJECT}_F_smooth250                 )
FIGURES=(  helpfiles/${PROJECT}/${PROJECT}_M.R    helpfiles/${PROJECT}/${PROJECT}_F.R    )

# ── Align reads (Step 2) ─────────────────────────────────────────────────────
NN=$(wc -l < ${BARCODES})
mkdir -p data/bam/${PROJECT}
jid_bam=$(sbatch --parsable --array=1-${NN} scripts/fq2bam.sh \
    ${BARCODES} data/raw/${PROJECT} data/bam/${PROJECT})
echo "fq2bam: $jid_bam"

# ── REFALT counts (Step 3) ───────────────────────────────────────────────────
mkdir -p process/${PROJECT}
find data/bam/${PROJECT} -name "*.bam" -size +1G > helpfiles/${PROJECT}/bams
cat helpfiles/founder.bams.txt >> helpfiles/${PROJECT}/bams

jid_refalt=$(sbatch --parsable --dependency=afterok:${jid_bam} \
    scripts/bam2bcf2REFALT.sh \
    helpfiles/${PROJECT}/bams process/${PROJECT})
echo "REFALT: $jid_refalt"

# ── Haplotypes (Step 4) ──────────────────────────────────────────────────────
jid_haps=$(sbatch --parsable --dependency=afterok:${jid_refalt} \
    --array=1-5 scripts/REFALT2haps.sh \
    --parfile ${PARFILE} --dir process/${PROJECT})
echo "haps:   $jid_haps"

# ── Haplotype scan + figures (Step 5a) ───────────────────────────────────────
for i in "${!DESIGNS[@]}"; do
    bash scripts/run_scan.sh \
        --design ${DESIGNS[$i]} \
        --dir    process/${PROJECT} \
        --scan   ${OUTDIRS[$i]} \
        --figure ${FIGURES[$i]} \
        --after  ${jid_haps}
done

# ── SNP scan (Step 5b, optional — remove this block if not needed) ───────────
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8
SNP_TABLE=helpfiles/FREQ_SNPs_Apop.cM.txt.gz

for i in "${!DESIGNS[@]}"; do
    bash scripts/run_snp_scan.sh \
        --design    ${DESIGNS[$i]} \
        --dir       process/${PROJECT} \
        --scan      ${OUTDIRS[$i]} \
        --snp-table ${SNP_TABLE} \
        --founders  ${FOUNDERS}
done
```

When all jobs finish, download results with:

```bash
scp <user>@<cluster>:<project_path>/process/<project>/<scan_name>/<scan_name>.tar.gz .
```

---

## Directory structure

```
XQTL2/
├── scripts/              # Core pipeline scripts (tracked in git)
├── scripts_oneoffs/      # Experiment-specific submit scripts (not tracked)
├── helpfiles/
│   ├── flymap.r6.txt
│   ├── founder.bams.txt
│   ├── FREQ_SNPs.cM.txt.gz              (SNP frequencies — see snp_tables/)
│   ├── FREQ_SNPs_Apop.cM.txt.gz         (A-pop subset, from prep_snp_table.R)
│   ├── FREQ_SNPs_Bpop.cM.txt.gz         (B-pop subset)
│   ├── snp_tables/README.md              (documents SNP table preparation)
│   └── <project>/
│       ├── <project>.barcodes.txt        (Step 2)
│       ├── bams                          (Step 3)
│       ├── hap_params.R                  (Step 4)
│       ├── design.txt                    (Step 5)
│       └── <figure_name>.R              (Step 6)
├── data/
│   ├── raw/<project>/                    (Step 1 — raw reads)
│   └── bam/<project>/                    (Step 2 — aligned bams)
├── ref/                  # Reference genome (not tracked)
├── process/
│   └── <project>/
│       ├── RefAlt.<chr>.txt              (Step 3)
│       ├── R.haps.<chr>.rds              (Step 4 — SNP table)
│       ├── R.haps.<chr>.out.rds          (Step 4 — haplotype estimates)
│       └── <scan_name>/                  (Step 5)
│           ├── <scan_name>.scan.txt
│           ├── <scan_name>.meansBySample.txt
│           ├── <scan_name>.snp_scan.txt
│           └── <scan_name>.tar.gz
└── figures/              # Publication figures (not tracked)
```

To check what exists for a given project on the cluster:

```bash
bash scripts/show_project_layout.sh <project>
```
