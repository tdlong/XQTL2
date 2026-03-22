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
   - 5a. Smooth haplotype frequencies (`smooth_haps.sh`)
   - 5b. Haplotype scan — Wald test + heritability (`hap_scan.sh`)
   - 5c. SNP scan — per-SNP Wald test (`snp_scan.sh`)
6. Concatenate chromosomes
7. Generate figures

A legacy scan without smoothing (`haps2scan.Apr2025.sh`) is also available as an
alternative to Steps 5a–5c.

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

## Step 5 — Run the scan

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

### 5a — Smooth haplotype frequencies

Applies a running-mean smoother (half-window = `--smooth-kb`) to haplotype
frequencies and covariance matrices. All downstream steps use these smoothed
estimates.

```bash
sbatch --array=1-5 scripts/smooth_haps.sh \
    --rfile     helpfiles/<project>/design.txt \
    --dir       process/<project> \
    --outdir    <scan_name> \
    --smooth-kb 125
```

### 5b — Haplotype scan (Wald test + heritability)

Runs a stabilized Wald test (eigenvalue-regularized) and two heritability
estimators (Falconer and Cutler) at each haplotype window.

```bash
sbatch --dependency=afterok:$JID_SMOOTH --array=1-5 scripts/hap_scan.sh \
    --rfile   helpfiles/<project>/design.txt \
    --dir     process/<project> \
    --outdir  <scan_name>
```

### 5c — SNP scan (optional, runs in parallel with 5b)

Imputes per-SNP allele frequencies from smoothed haplotype estimates, then runs a
Wald test (df=1) at every SNP. Requires a pre-built SNP table — see
`helpfiles/snp_tables/` for how to prepare one.

```bash
sbatch --dependency=afterok:$JID_SMOOTH --array=1-5 scripts/snp_scan.sh \
    --rfile     helpfiles/<project>/design.txt \
    --dir       process/<project> \
    --outdir    <scan_name> \
    --snp-table helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founders  A1,A2,A3,A4,A5,A6,A7,AB8
```

### Per-chromosome outputs

Each scan step writes one file per chromosome:

```
process/<project>/<scan_name>/
    <scan_name>.scan.<chr>.txt              (from hap_scan)
    <scan_name>.meansBySample.<chr>.txt     (from smooth_haps)
    <scan_name>.snp_scan.<chr>.txt          (from snp_scan)
```

### Legacy scan (alternative — no smoothing)

```bash
sbatch --array=1-5 scripts/haps2scan.Apr2025.sh \
    --rfile  helpfiles/<project>/design.txt \
    --dir    process/<project> \
    --outdir <scan_name>
```

---

## Step 6 — Concatenate chromosomes

### Haplotype scan

```bash
bash scripts/concat_Chromosome_Scans.sh process/<project>/<scan_name>
```

Merges per-chromosome scan and meansBySample files, generates quick-look Manhattan plots,
and bundles everything into `<scan_name>.tar.gz`.

### SNP scan

```bash
bash scripts/concat_snp_scans.sh process/<project>/<scan_name>
```

Merges per-chromosome `snp_scan` files into `<scan_name>.snp_scan.txt`.

---

## Step 7 — Generate publication figures

Three plot engines are available, each driven by a small R script that sets parameters
then `source()`s the engine. Save figure scripts to `helpfiles/<project>/`.

### Wald scan — `plot_pseudoscan.R`

5-panel -log10(p) line plot from the haplotype scan.

```r
SCAN_FILES   <- c("process/<project>/<scan_name>/<scan_name>.scan.txt")
SCAN_COLOURS <- c("#1F78B4")
SCAN_LABELS  <- NULL       # NULL uses the file basename
THRESHOLD    <- 10         # dashed horizontal line
OUT_FILE     <- "process/<project>/<scan_name>/<figure_name>.png"
FORMAT       <- "powerpoint"
PEAKS        <- NULL
GENES        <- NULL

source("scripts/plot_pseudoscan.R")
```

Two-scan overlay (e.g. males and females):

```r
SCAN_FILES   <- c("path/to/<scan_M>.scan.txt", "path/to/<scan_F>.scan.txt")
SCAN_COLOURS <- c("#1F78B4", "#E31A1C")
SCAN_LABELS  <- c("Male", "Female")
THRESHOLD    <- 10
OUT_FILE     <- "path/to/overlay.png"
FORMAT       <- "powerpoint"
PEAKS        <- NULL
GENES        <- NULL

source("scripts/plot_pseudoscan.R")
```

### Heritability overlay — `plot_H2_overlay.R`

5-panel line plot overlaying Falconer H2 and Cutler H2 from the haplotype scan.

```r
SCAN_FILE <- "process/<project>/<scan_name>/<scan_name>.scan.txt"
OUT_FILE  <- "process/<project>/<scan_name>/<figure_name>.H2.png"
FORMAT    <- "powerpoint"

source("scripts/plot_H2_overlay.R")
```

### SNP scan — `plot_freqsmooth_snp.R`

5-panel -log10(p) dot plot from the SNP scan.

```r
SCAN_FILE  <- "process/<project>/<scan_name>/<scan_name>.snp_scan.txt"
OUT_FILE   <- "process/<project>/<scan_name>/<figure_name>.snp.png"
FORMAT     <- "powerpoint"
THRESHOLD  <- 10

source("scripts/plot_freqsmooth_snp.R")
```

### Run a figure script

```bash
module load R/4.2.2
Rscript helpfiles/<project>/<figure_name>.R
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

Gene and peak annotations use Mb coordinates:

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

## Step 8 — Download and explore results

The haplotype concat (Step 6) bundles scan tables and quick-look figures into a tarball:

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

## Putting it all together

For a new experiment, copy the block below to `scripts_oneoffs/<project>_pipeline.sh`,
fill in the variables, and run it. Each step is submitted as a SLURM job with
`afterok` dependencies so the whole pipeline runs unattended.

```bash
#!/bin/bash
# Edit these variables, then: bash scripts_oneoffs/<project>_pipeline.sh

PROJECT=<project>
BARCODES=helpfiles/${PROJECT}/${PROJECT}.barcodes.txt
PARFILE=helpfiles/${PROJECT}/hap_params.R
SMOOTH_KB=125
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8
SNP_TABLE=helpfiles/FREQ_SNPs_Apop.cM.txt.gz

# One entry per scan (sex, treatment, etc.)
DESIGNS=(  helpfiles/${PROJECT}/design_male.txt   helpfiles/${PROJECT}/design_female.txt )
OUTDIRS=(  ${PROJECT}_M_smooth${SMOOTH_KB}        ${PROJECT}_F_smooth${SMOOTH_KB}        )

# Figure scripts — run after all scans are concatenated
FIGURES=(  helpfiles/${PROJECT}/${PROJECT}_MF.R )

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
    scripts/bam2bcf2REFALT.sh \
    helpfiles/${PROJECT}/bams process/${PROJECT})
echo "REFALT:      $jid_refalt"

# ── Step 4: haplotypes ────────────────────────────────────────────────────────
jid_haps=$(sbatch --parsable --dependency=afterok:${jid_refalt} \
    --array=1-5 scripts/REFALT2haps.sh \
    --parfile ${PARFILE} --dir process/${PROJECT})
echo "haps:        $jid_haps"

# ── Steps 5–7: scan + concat + figures (one set per design) ──────────────────
jid_concats=""
for i in "${!DESIGNS[@]}"; do
    design=${DESIGNS[$i]}
    outdir=${OUTDIRS[$i]}

    # 5a: smooth
    jid_smooth=$(sbatch --parsable --dependency=afterok:${jid_haps} \
        --array=1-5 scripts/smooth_haps.sh \
        --rfile ${design} --dir process/${PROJECT} \
        --outdir ${outdir} --smooth-kb ${SMOOTH_KB})
    echo "smooth ${outdir}: $jid_smooth"

    # 5b: haplotype scan
    jid_scan=$(sbatch --parsable --dependency=afterok:${jid_smooth} \
        --array=1-5 scripts/hap_scan.sh \
        --rfile ${design} --dir process/${PROJECT} --outdir ${outdir})
    echo "hap_scan ${outdir}: $jid_scan"

    # 5c: SNP scan (parallel with 5b)
    jid_snp=$(sbatch --parsable --dependency=afterok:${jid_smooth} \
        --array=1-5 scripts/snp_scan.sh \
        --rfile ${design} --dir process/${PROJECT} --outdir ${outdir} \
        --snp-table ${SNP_TABLE} --founders ${FOUNDERS})
    echo "snp_scan ${outdir}: $jid_snp"

    # 6: concatenate
    jid_concat=$(sbatch --parsable --dependency=afterok:${jid_scan} \
        -A tdlong_lab -p standard --mem=10G \
        --wrap="bash scripts/concat_Chromosome_Scans.sh process/${PROJECT}/${outdir}")
    echo "concat ${outdir}: $jid_concat"

    jid_snp_concat=$(sbatch --parsable --dependency=afterok:${jid_snp} \
        -A tdlong_lab -p standard --mem=10G \
        --wrap="bash scripts/concat_snp_scans.sh process/${PROJECT}/${outdir}")
    echo "snp_concat ${outdir}: $jid_snp_concat"

    jid_concats="${jid_concats}:${jid_concat}:${jid_snp_concat}"
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
│       └── <figure_name>.R              (Step 7)
├── data/
│   ├── raw/<project>/                    (Step 1 — raw reads)
│   └── bam/<project>/                    (Step 2 — aligned bams)
├── ref/                  # Reference genome (not tracked)
├── process/
│   └── <project>/
│       ├── RefAlt.<chr>.txt              (Step 3)
│       ├── R.haps.<chr>.rds              (Step 4 — SNP table)
│       ├── R.haps.<chr>.out.rds          (Step 4 — haplotype estimates)
│       └── <scan_name>/                  (Steps 5–6)
│           ├── <scan_name>.scan.<chr>.txt
│           ├── <scan_name>.meansBySample.<chr>.txt
│           ├── <scan_name>.snp_scan.<chr>.txt
│           ├── <scan_name>.scan.txt          (after concat)
│           ├── <scan_name>.meansBySample.txt
│           ├── <scan_name>.snp_scan.txt
│           └── <scan_name>.tar.gz
└── figures/              # Publication figures (not tracked)
```

To check what exists for a given project on the cluster:

```bash
bash scripts/show_project_layout.sh <project>
```
