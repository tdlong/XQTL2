# XQTL Pipeline

## Overview

This pipeline takes pooled-sequencing XQTL data from raw reads to genome scans
and publication figures. Everything runs on a SLURM cluster. You can submit the
entire pipeline with a single script and walk away, or run each step individually
for greater control.

Once you have haplotype calls you can run additional scans with different designs
and generate new figures without repeating the expensive alignment and haplotype
steps. We are also developing [XQTL.xplore](https://github.com/tdlong/XQTL.xplore),
a companion package for interactive graphical analysis of scan results.

**Pipeline at a glance:**

1. Get raw reads
2. Align reads (`fq2bam.sh`)
3. Generate allele counts (`bam2bcf2REFALT.sh`)
4. Call haplotypes (`REFALT2haps.sh`)
5. Scan
   - 5a. Haplotype scan (`run_scan.sh` — smooth, Wald test, H², concat)
   - 5b. SNP scan (`run_snp_scan.sh` — optional, imputed SNP-level test)
6. Generate figures (`plot_5panel.R`, `plot_manhattan.R`, `plot_H2_overlay.R`, `plot_freqsmooth_snp.R`)
7. Download results

**What's already in this repo:**

- All pipeline scripts (`scripts/`)
- Founder BAM paths (`helpfiles/A_founders.bams.txt`) — paths relative to this repo
- Per-founder SNP state tables for SNP frequency imputation (`helpfiles/FREQ_SNPs_Apop.cM.txt.gz`, `FREQ_SNPs_Bpop.cM.txt.gz`)
- Physical-to-genetic map (`helpfiles/flymap.r6.txt`)
- Heterochromatic boundary definitions (`helpfiles/het_bounds.txt`)
- Generic haplotype parameter template (`helpfiles/generic_haplotype_parameters.R`)

**What you need to provide per experiment:**

- Raw sequencing reads (from the sequencing core)
- A barcode file mapping barcodes → sample names (Step 2)
- A haplotype parameters file listing your founders, samples, and window sizes (Step 4)
- A design file describing your experimental layout (Step 5)

**What you get at the end:**

- `*.scan.txt` — genome-wide haplotype scan (Wald -log10p, Falconer H², Cutler H²)
- `*.snp_scan.txt` — SNP-level Wald test at every imputed SNP (optional)
- `*.meansBySample.txt` — smoothed founder frequencies per sample (QC)
- Manhattan plots and heritability figures (PNG)
- Tarballs of everything, ready to scp down

The [mylab-XQTL template](https://github.com/tdlong/mylab-XQTL) contains a
fully configured malathion resistance experiment used as a training example — see
the **Worked example** section at the end of this README.

---

## Installation and Setup

This is a one-time setup per machine. You need two things: the pipeline (this
repo) and a project repo for your own data and scripts.

### 1. Clone the pipeline

```bash
git clone https://github.com/tdlong/XQTL2.git
cd XQTL2
```

### 2. Download founder BAMs

Pre-aligned founder BAMs are hosted on the Long lab server. Download and unpack
into `data/founders/`:

```bash
mkdir -p data/founders
wget https://wfitch.bio.uci.edu/~tdlong/founders_bam_files.tar
tar -xf founders_bam_files.tar -C data/founders/
rm founders_bam_files.tar
```

The founder BAM paths are pre-configured in `helpfiles/A_founders.bams.txt` and `helpfiles/B_founders.bams.txt` — no editing needed.

### 3. Download and index the reference genome

The pipeline aligns to dm6 (*Drosophila melanogaster* release 6). We use
UCSC-style chromosome names (`chr2L`, `chr3R`, etc.) rather than FlyBase names.
Download from UCSC:

```bash
mkdir -p ref
wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip dm6.fa.gz
mv dm6.fa ref/
bwa index ref/dm6.fa
samtools faidx ref/dm6.fa
samtools dict ref/dm6.fa > ref/dm6.dict
```

Submit as a SLURM job — BWA indexing takes ~1 hour. The pipeline is not limited
to Drosophila; it can be adapted to any synthetic population with pooled-sequencing XQTL data
by substituting the appropriate reference genome, genetic map, and founder BAMs.

### 4. Create your project repo

Your experimental data and project-specific scripts live in a separate repo
alongside XQTL2. The [mylab-XQTL template](https://github.com/tdlong/mylab-XQTL)
gives you a ready-made starting point — clone it, rename it, and make it your
own. This is one of the key advantages of separating your project from the
pipeline: you get your own git repo from day one.

```bash
cd ..   # move up next to XQTL2/
git clone https://github.com/tdlong/mylab-XQTL.git LongLab-XQTL
cd LongLab-XQTL
ln -s ../XQTL2 pipeline   # all pipeline calls go through this symlink
```

Replace `LongLab-XQTL` with whatever makes sense for your group. The
`pipeline` symlink is the only path that matters.

All your submission scripts call `pipeline/scripts/run_scan.sh` etc. through
the symlink. When XQTL2 is updated, run `git pull` inside the XQTL2 directory
— your scripts automatically use the new version.

### Version control with git (optional but recommended)

Git is not required — the pipeline works fine without it. But because your
project is now its own directory, separate from the pipeline, it is easy to
put it under version control. If you do, you get a complete record of every
config file you created and every change you made, and you can sync your
project across machines (e.g. keep your laptop in sync with the cluster).

**One-time git setup** (skip if you've used git before):

```bash
git config --global user.name "Your Name"
git config --global user.email "you@example.com"
```

You also need a way to authenticate with GitHub. The simplest option is to
install the [GitHub CLI](https://cli.github.com/) and run `gh auth login`,
which walks you through it interactively. Alternatively, you can create a
[personal access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens)
on GitHub and paste it when prompted for a password.

**Create your GitHub repo:** go to [github.com/new](https://github.com/new),
name it whatever you named your directory (e.g. `LongLab-XQTL`), and leave
it empty (no README, no .gitignore — your clone already has these). Then
point your local clone at it and push:

```bash
git remote set-url origin https://github.com/<you>/<YourLab-XQTL>.git
git push -u origin main
```

From here on, the git workflow is just four commands:

```bash
git add <files>                       # stage new or changed files
git commit -m "describe what changed" # save a snapshot
git push                              # upload to GitHub
git pull                              # download changes on another machine
```

The template's `.gitignore` already excludes large files (BAMs, scan output,
raw reads). What gets tracked is everything *you* created: barcode files, BAM
lists, haplotype parameters, design files, and submission scripts. As long as
the raw FASTQs are backed up separately, every result can be regenerated from
these tracked config files plus the pipeline.

The steps below will note where it makes sense to commit.

---

### SLURM resource requirements

The cluster's standard partition provides max 6 GB per core; highmem provides
10 GB per core. Always use `--mem-per-cpu` (not `--mem`). See `Slurm.md` for
full partition details.

Default SLURM account and partition are `tdlong_lab` and `standard`. Pass
`-A <account>` and `-p <partition>` to any pipeline script to override — all
scripts accept these flags.

Scan steps were profiled with `seff` on the malathion training dataset (2
effective replicates, 4 samples). Larger experiments scale proportionally.

| Script | Step | Partition | CPUs | Mem/CPU | Time | Profiled (malathion) |
|--------|------|-----------|------|---------|------|----------------------|
| `fq2bam.sh` | 2 | standard | 4 | 6G | 1 day | `bwa -t 4`; `java -Xmx20g` needs ~20G total |
| `bam2bcf2REFALT.sh` | 3 | standard | 2 | 6G | 5 days | bcftools mpileup, I/O-bound |
| `REFALT2haps.sh` | 4 | highmem | 1 | 10G | 1 day | large haplotype matrices require highmem |
| `smooth_haps.sh` | 5a | standard | 1 | 3G | 4 hr | 909 MB / 17s wall |
| `hap_scan.sh` | 5a | standard | 1 | 3G | 4 hr | 307 MB / 5:12 wall |
| `snp_scan.sh` | 5b | standard | 1 | 3G | 4 hr | 732 MB / 5:25 wall |
| concat | 5a | standard | 1 | 3G | 1 hr | 436 MB / 19s wall |
| snp_concat | 5b | standard | 1 | 3G | 1 hr | 413 MB / 21s wall |
| figures | 6 | standard | 1 | 3G | 1 hr | 982 MB / 58s wall |

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

### Barcode-to-sample mapping file

Create a tab-delimited file mapping sequencing barcodes to sample names. Each row
is one sample: forward barcode, reverse barcode, sample name. Sample names become
the BAM file prefixes and read group IDs used throughout the pipeline.

```
TGGCTATG    TTGTCAGC    R3con
GTCCTAGA    TTGTCAGC    R3age
ACTTGCCA    TTGTCAGC    R5con
TCTTCGTG    TTGTCAGC    R5age
```

Save this file to `helpfiles/<project>/<project>.barcodes.txt`. This is your
first project-specific config file — a good time to commit:

```bash
git add helpfiles/<project>/<project>.barcodes.txt
git commit -m "add barcode file for <project>"
git push
```

### Run alignment

```bash
mkdir -p data/bam/<project>
NN=$(wc -l < helpfiles/<project>/<project>.barcodes.txt)
sbatch --array=1-$NN pipeline/scripts/fq2bam.sh \
    helpfiles/<project>/<project>.barcodes.txt \
    data/raw/<project> \
    data/bam/<project>
```

`fq2bam.sh` handles barcode splitting, BWA alignment, coordinate sorting, and
Picard AddOrReplaceReadGroups in one job per sample. **Sample names in the output
BAMs are set by Picard read groups and must exactly match what downstream scripts
expect.** Verify with `samtools view -H <sample>.bam | grep "^@RG"`.

Bam files below ~1 GB likely indicate a failed library prep and should be reprocessed.

---

## Step 3 — Generate REFALT counts (bam to REFALT)

Create `helpfiles/<project>/bam_list.txt` — one BAM path per line, sample BAMs
first, then founders. Build a draft from your BAM directory, append founders,
then **review it** before submitting.

First, check which founders are available:

```bash
cat pipeline/helpfiles/A_founders.bams.txt
# pipeline/data/founders/A1.dedup.bam
# pipeline/data/founders/A2.dedup.bam
# pipeline/data/founders/A3.dedup.bam
# pipeline/data/founders/A4.dedup.bam
# pipeline/data/founders/A5.dedup.bam
# pipeline/data/founders/A6.dedup.bam
# pipeline/data/founders/A7.dedup.bam
# pipeline/data/founders/AB8.dedup.bam

cat pipeline/helpfiles/B_founders.bams.txt
# pipeline/data/founders/AB8.dedup.bam
# pipeline/data/founders/B1.dedup.bam
# ...
# pipeline/data/founders/B7.dedup.bam
```

Use the founder file that matches your population. If your experiment used A-pop
founders, use `A_founders.bams.txt` (A1–A7, AB8). If B-pop, use `B_founders.bams.txt`
(AB8, B1–B7). Your `hap_params.R` founders list must match exactly.

If your design crossed the synthetic population to a tester strain or other
reference genotype, treat that strain as an additional founder and include its
BAM in the list alongside the population founders.

```bash
# Draft from your sample BAMs
ls data/bam/<project>/*.bam > helpfiles/<project>/bam_list.txt

# Append the founders that match your experiment
cat pipeline/helpfiles/A_founders.bams.txt >> helpfiles/<project>/bam_list.txt

# Review — confirm every sample and every founder is present, no extras
cat helpfiles/<project>/bam_list.txt

# Commit and push — this records exactly what went into your analysis
git add helpfiles/<project>/bam_list.txt
git commit -m "add bam list for <project>"
git push
```

```bash
# Submit
mkdir -p process/<project>
sbatch --array=1-5 pipeline/scripts/bam2bcf2REFALT.sh \
    helpfiles/<project>/bam_list.txt \
    process/<project>
```

This runs as a 5-task array (one chromosome each: chrX, chr2L, chr2R, chr3L,
chr3R) and produces `RefAlt.<chr>.txt` files in `process/<project>/`.

---

## Step 4 — Call haplotypes (REFALT to haps)

### Haplotype parameters file

Create `helpfiles/<project>/hap_params.R`. Founder names must exactly match the
read group sample names in the founder BAMs (set by Picard during alignment):

```r
# Founder set — names must match read groups in the founder BAMs
founders <- c("A1","A2","A3","A4","A5","A6","A7","AB8")

# Sample names — must exactly match BAM prefixes from Step 2
names_in_bam <- c("R1con","R1age","R2con","R2age","R3con","R3age",
                   "R4con","R4age","R5con","R5age","R6con","R6age")

# Window step size in bp (5000 typical; 10000 for very large experiments)
step <- 5000

# Base half-window in bp. In low-recombination regions the window expands
# proportional to max_RR / local_RR so each window spans similar
# recombination distances regardless of genomic position.
size <- 50000

# Tree height cutoff for founder distinguishability (2.5 is default)
h_cutoff <- 2.5
```

See `pipeline/helpfiles/generic_haplotype_parameters.R` for a full template with
comments. To generate `names_in_bam` from your BAM directory:

```bash
echo -n "names_in_bam <- c(" && \
find data/bam/<project> -name "*.bam" -size +1G -print0 | \
xargs -0 -n1 basename | sed 's/.bam//' | sort | \
sed 's/.*/"&"/' | tr '\n' ',' | sed 's/,$//' && echo ")"
```

Commit the parameters file:

```bash
git add helpfiles/<project>/hap_params.R
git commit -m "add haplotype parameters for <project>"
git push
```

### Run haplotype calling

This runs as a 5-task array, one chromosome per task:

```bash
sbatch --array=1-5 pipeline/scripts/REFALT2haps.sh \
    --parfile helpfiles/<project>/hap_params.R \
    --dir     process/<project>
```

---

## Step 5 — Scan

The scan takes the haplotype output and runs the Wald test, heritability
estimates, and chromosome concatenation. There are two ways to run it:

- **`run_scan.sh` (recommended):** submits smoothing → scan → concat as a
  chained SLURM pipeline with one command.
- **Step by step:** submit each stage individually — useful if you want to
  adjust parameters between stages or diagnose failures.

### Design file

Create a plain text table with one row per sample. **Column names are
case-sensitive** — the pipeline refers to them by exact name.

| Column | Description |
|--------|-------------|
| `bam` | Sample name (must match BAM prefix from Step 2) |
| `TRT` | `C` = control, `Z` = selected |
| `REP` | Replicate number (integer) |
| `REPrep` | Technical replicate within replicate (usually `1`) |
| `Num` | Number of flies in pool |
| `Proportion` | Fraction selected (`NA` for controls) |

Create and save from R (use `write.table` defaults — row numbers are included):

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

```bash
git add helpfiles/<project>/design.txt
git commit -m "add design file for <project>"
git push
```

### Option A — One-command scan (recommended)

`run_scan.sh` chains smoothing → Wald test + H² → chromosome concat with proper
SLURM dependency chaining. One command per scan:

```bash
bash pipeline/scripts/run_scan.sh \
    --design    helpfiles/<project>/design.txt \
    --dir       process/<project> \
    --scan      <scan_name> \
    --after     $JID_HAPS
```

| Flag | Default | Description |
|------|---------|-------------|
| `--design` | (required) | Path to design file |
| `--dir` | (required) | Project directory (e.g. `process/<project>`) |
| `--scan` | (required) | Scan name — becomes output subdirectory |
| `--smooth` | 250 | Smoothing half-window in kb |
| `--mem-per-cpu` | 3G | Memory per CPU for all jobs |
| `--cpus-per-task` | 1 | CPUs for all jobs |
| `-p` / `--partition` | standard | SLURM partition |
| `-A` / `--account` | tdlong_lab | SLURM account |
| `--after` | (none) | SLURM job ID to wait on (e.g. from REFALT2haps) |

### Option B — Step by step

For reference or debugging, `run_scan.sh` internally chains these three jobs:

**1. Smooth haplotype frequencies** (5-task array, one chromosome per task):

```bash
sbatch --array=1-5 pipeline/scripts/smooth_haps.sh \
    --rfile     helpfiles/<project>/design.txt \
    --dir       process/<project> \
    --outdir    <scan_name> \
    --smooth-kb 250
```

**2. Haplotype scan** — Wald test + heritability at each window (5-task array,
after smooth step):

```bash
sbatch --array=1-5 pipeline/scripts/hap_scan.sh \
    --rfile  helpfiles/<project>/design.txt \
    --dir    process/<project> \
    --outdir <scan_name>
```

**3. Concatenate chromosomes and generate quick-look plots:**

```bash
bash pipeline/scripts/concat_scans.sh process/<project>/<scan_name>
```

### Smoothing window

The default smoothing window is 250 kb, chosen by comparing 125 kb and 250 kb
on real data (malathion experiment). At 250 kb the founder frequency estimates
are stable without over-smoothing genuine biological signal. Override with
`--smooth 125` or any value.

### Legacy scan (no smoothing)

An older scan method without haplotype smoothing is preserved for reference:

```bash
sbatch --array=1-5 pipeline/scripts/haps2scan.Apr2025.sh \
    --rfile  helpfiles/<project>/design.txt \
    --dir    process/<project> \
    --outdir <scan_name>
bash pipeline/scripts/concat_scans.sh process/<project>/<scan_name>
```

---

## Step 5b — SNP scan (optional)

The SNP scan imputes per-SNP ALT allele frequencies from the smoothed haplotype
estimates and runs a Wald test (df=1) at every SNP. It tests at individual SNP
positions rather than haplotype windows, but the signal comes from the same
smoothed haplotype estimates — it is not independent of the haplotype scan.

`run_scan.sh` must have already run before starting the SNP scan.

### SNP table

The scan requires per-founder SNP state tables — one-time preparation per
population. See `pipeline/helpfiles/snp_tables/` for preparation details.
Pre-built tables for A-pop and B-pop are included in this repo.

### Run the SNP scan

```bash
bash pipeline/scripts/run_snp_scan.sh \
    --design    helpfiles/<project>/design.txt \
    --dir       process/<project> \
    --scan      <scan_name> \
    --snp-table pipeline/helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founders  A1,A2,A3,A4,A5,A6,A7,AB8
```

Use the same `--scan` name as Step 5a. Options mirror `run_scan.sh`.

---

## Step 6 — Generate publication figures

Three plotting scripts produce 5-panel (per-chromosome) figures. These are
submitted via `--wrap` in the worked example scripts, or run interactively.

### Haplotype Wald scan

```bash
Rscript pipeline/scripts/plot_5panel.R \
    --scan   process/<project>/<scan_name>/<scan_name>.scan.txt \
    --out    process/<project>/<scan_name>/<scan_name>.wald.png \
    --format powerpoint \
    --threshold 10
```

Overlay two scans (e.g. male vs female):

```bash
Rscript pipeline/scripts/plot_5panel.R \
    --scan   process/<project>/<scan_M>/<scan_M>.scan.txt \
    --scan   process/<project>/<scan_F>/<scan_F>.scan.txt \
    --label  Male --label Female \
    --colour "#1F78B4" --colour "#E31A1C" \
    --out    process/<project>/MF_overlay.png \
    --format powerpoint --threshold 10
```

### Manhattan plot (single-row)

A traditional Manhattan with all chromosomes concatenated on one x-axis.
Heterochromatic regions are shaded and chromosome boundaries marked with
dotted vertical lines. Same interface as `plot_5panel.R`:

```bash
Rscript pipeline/scripts/plot_manhattan.R \
    --scan   process/<project>/<scan_name>/<scan_name>.scan.txt \
    --out    process/<project>/<scan_name>/<scan_name>.manhattan.png \
    --format powerpoint \
    --threshold 10
```

### Heritability overlay (Falconer + Cutler)

```bash
Rscript pipeline/scripts/plot_H2_overlay.R \
    --scan   process/<project>/<scan_name>/<scan_name>.scan.txt \
    --out    process/<project>/<scan_name>/<scan_name>.H2.png \
    --format powerpoint
```

### SNP scan

```bash
Rscript pipeline/scripts/plot_freqsmooth_snp.R \
    --scan   process/<project>/<scan_name>/<scan_name>.snp_scan.txt \
    --out    process/<project>/<scan_name>/<scan_name>.snp.wald.png \
    --format powerpoint --threshold 10
```

### Common options

| Flag | Description |
|------|-------------|
| `--scan <file>` | Input scan file (required; repeat for overlays) |
| `--out <file>` | Output PNG path (required) |
| `--format <name>` | Size/DPI preset (default: `powerpoint`) |
| `--threshold <n>` | Dashed horizontal line at this y value |
| `--genes <file>` | Tab-delimited gene annotations (`name`, `chr`, `pos_mb`) |
| `--peaks <file>` | Tab-delimited peak annotations (`label`, `chr`, `pos_mb`) |
| `--height <in>` | Override figure height in inches |

`plot_5panel.R` and `plot_freqsmooth_snp.R` also accept `--label` and
`--colour` (one per `--scan`) for overlays.

### FORMAT presets

| FORMAT | Width | DPI | Use for |
|--------|-------|-----|---------|
| `manuscript_half` | 3.5 in | 300 | half-width journal figure |
| `manuscript_full` | 7.0 in | 300 | full-width journal figure |
| `manuscript_half_hires` | 3.5 in | 600 | high-res submission |
| `manuscript_full_hires` | 7.0 in | 600 | high-res submission |
| `powerpoint` | 8.0 in | 150 | slides |
| `web` | 7.0 in | 150 | web/HTML |
| `email` | 6.0 in | 100 | email preview |

### Gene and peak annotations

```
# helpfiles/<project>/genes.txt
name	chr	pos_mb
Ace	chr3R	9.07
Cyp6g1	chr2R	12.19
```

Pass with `--genes helpfiles/<project>/genes.txt` (works with all three plotters).

### Interactive exploration (optional)

```r
source("pipeline/scripts/XQTL_plotting_functions.R")
df1 <- as_tibble(read.table("process/<project>/<scan>/<scan>.scan.txt"))
df2 <- as_tibble(read.table("process/<project>/<scan>/<scan>.meansBySample.txt"))

XQTL_Manhattan_5panel(df1, cM = FALSE)
XQTL_region(df1, "chr3R", 18250000, 19000000, "Wald_log10p")
XQTL_change_average(df2, "chr3R", 18250000, 19000000)
XQTL_combined_plot(df1, df2, "chr3R", 18250000, 19000000)
```

---

## End-to-end: from raw reads to figures

This section walks through a complete experiment from start to finish. If you
are using an AI assistant, point it at this README — it can generate the
commands below for your specific project.

Suppose your experiment is called `heatshock`, uses A-population founders, and
has 6 samples (3 control, 3 selected). You've already completed Installation
steps 1–4 and downloaded your raw reads to `data/raw/heatshock/`.

### 1. Create your config files

These are the four files the pipeline needs. The steps above (Steps 2–5)
explain each one in detail — here we just show the end result.

```bash
mkdir -p helpfiles/heatshock
```

**Barcode file** (`helpfiles/heatshock/heatshock.barcodes.txt`) — maps
sequencing barcodes to sample names:

```
TGGCTATG	TTGTCAGC	R1con
GTCCTAGA	TTGTCAGC	R1heat
ACTTGCCA	TTGTCAGC	R2con
TCTTCGTG	TTGTCAGC	R2heat
AAGCGACT	TTGTCAGC	R3con
CGTGAATC	TTGTCAGC	R3heat
```

**Haplotype parameters** (`helpfiles/heatshock/hap_params.R`):

```r
founders     <- c("A1","A2","A3","A4","A5","A6","A7","AB8")
names_in_bam <- c("R1con","R1heat","R2con","R2heat","R3con","R3heat")
step         <- 5000
size         <- 50000
h_cutoff     <- 2.5
```

**Design file** (`helpfiles/heatshock/design.txt`) — create in R:

```r
design <- data.frame(
    bam        = c("R1con","R1heat","R2con","R2heat","R3con","R3heat"),
    TRT        = c("C","Z","C","Z","C","Z"),
    REP        = c(1,1,2,2,3,3),
    REPrep     = 1,
    Num        = c(500,100,500,100,500,100),
    Proportion = c(NA,0.20,NA,0.20,NA,0.20)
)
write.table(design, "helpfiles/heatshock/design.txt")
```

**Commit everything before running:**

```bash
git add helpfiles/heatshock/
git commit -m "add config files for heatshock experiment"
git push
```

### 2. Run the full pipeline

`run_full_pipeline.sh` chains every step (align → REFALT → haplotypes →
scan → SNP scan → figures) with SLURM dependency chaining. Submit once and
walk away:

```bash
bash pipeline/scripts/run_full_pipeline.sh \
    --project      heatshock \
    --barcodes     helpfiles/heatshock/heatshock.barcodes.txt \
    --rawdir       data/raw/heatshock \
    --bamdir       data/bam/heatshock \
    --parfile      helpfiles/heatshock/hap_params.R \
    --design       helpfiles/heatshock/design.txt \
    --scan         heatshock_smooth250 \
    --founders     A \
    --snp-table    pipeline/helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founder-list A1,A2,A3,A4,A5,A6,A7,AB8
```

The script prints each SLURM job ID as it submits. When the final job
finishes, results are in `process/heatshock/heatshock_smooth250/`.

### 3. Download results

```bash
scp <user>@<cluster>:<path>/process/heatshock/heatshock_smooth250/heatshock_smooth250.tar.gz .
tar xzf heatshock_smooth250.tar.gz
```

### Rerunning with different parameters

If you already have BAMs and want to rerun from REFALT (e.g. different
window sizes):

```bash
bash pipeline/scripts/run_full_pipeline.sh \
    --skip-fq2bam \
    --project heatshock --parfile helpfiles/heatshock/hap_params.R \
    --design helpfiles/heatshock/design.txt --scan heatshock_v2 \
    --founders A --snp-table pipeline/helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founder-list A1,A2,A3,A4,A5,A6,A7,AB8
```

If you already have haplotypes and just want a new scan with a different
design (e.g. different contrasts or subsets of samples):

```bash
bash pipeline/scripts/run_full_pipeline.sh \
    --skip-fq2bam --skip-refalt \
    --project heatshock --parfile helpfiles/heatshock/hap_params.R \
    --design helpfiles/heatshock/design_subset.txt --scan heatshock_subset \
    --founders A --snp-table pipeline/helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founder-list A1,A2,A3,A4,A5,A6,A7,AB8
```

All SLURM flags (`--mem-per-cpu`, `-p`, `-A`) are passed through to every job.

---

## Step 7 — Download results

The figure scripts bundle everything into a tarball at the end of the worked
example pipeline. Download it:

```bash
scp <user>@<cluster>:<path>/process/<project>/<scan_name>/<scan_name>.tar.gz .
tar xzf <scan_name>.tar.gz
```

**Tarball contents:**

| File | Contents |
|------|----------|
| `<scan>.scan.txt` | Haplotype scan (Wald -log10p, H²; one row per window) |
| `<scan>.meansBySample.txt` | Smoothed founder frequencies per window × treatment × rep × founder |
| `<scan>.snp_scan.txt` | SNP scan (Wald -log10p; one row per SNP) |
| `<scan>.snp_meansBySample.txt` | Imputed SNP ALT frequencies |
| `<scan>.wald.png` | 5-panel haplotype Wald Manhattan |
| `<scan>.H2.png` | 5-panel Falconer + Cutler heritability overlay |
| `<scan>.snp.wald.png` | 5-panel SNP Wald Manhattan |

**Output column reference:**

`<scan>.scan.txt` (one row per haplotype window):
`chr`, `pos` (bp), `Wald_log10p`, `Falc_H2`, `Cutl_H2`, `cM`

`<scan>.snp_scan.txt` (one row per SNP):
`chr`, `pos`, `Wald_log10p`, `cM`, `n_informative_founders`

`<scan>.meansBySample.txt`:
`chr`, `pos`, `TRT`, `REP`, `founder`, `freq`

---

## Worked example — malathion resistance

The [mylab-XQTL template](https://github.com/tdlong/mylab-XQTL) includes a
complete malathion resistance experiment pre-configured as a training run.
Running it end to end verifies that your pipeline installation is working
correctly.

**Experimental design:** Four pools of *D. melanogaster* (A-population) were
exposed to malathion and the survivors sequenced alongside untreated controls —
two sexes, one replicate each. The dataset is intentionally small (43–65 flies
per pool) so the full pipeline finishes quickly.

**Setup** (after completing Installation steps 1–4):

```bash
# Download the malathion BAMs into your project repo
mkdir -p data/bam/malathion
wget https://wfitch.bio.uci.edu/~tdlong/malathion_bams.tar
tar -xf malathion_bams.tar -C data/bam/malathion/
rm malathion_bams.tar
```

**Run the pipeline** (starts at Step 3 — BAMs are already provided):

```bash
bash pipeline/scripts/run_full_pipeline.sh \
    --skip-fq2bam \
    --project     malathion \
    --parfile     helpfiles/malathion/hap_params.R \
    --design      helpfiles/malathion/design.txt \
    --scan        malathion_smooth250 \
    --founders    A \
    --snp-table   pipeline/helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founder-list A1,A2,A3,A4,A5,A6,A7,AB8
```

This submits Steps 3–6 with SLURM dependency chaining and prints each job ID.
When complete, download the results tarball and check for a signal on chr3R
around 9 Mb. For biological interpretation see [Long et al. 2022](https://pubmed.ncbi.nlm.nih.gov/36250804/).

This is exactly the kind of command an AI assistant can generate for your own
experiment — point it at this README and your config files, and ask it to write
the `run_full_pipeline.sh` invocation for your project.

---

## Worked example — adding replicates to an existing experiment

You sequenced 3 replicates, ran the pipeline, then sequenced 3 more. You only
need to align the new samples, then rerun from Step 3 with all BAMs combined.

**What changes:**
- New barcode file for the new samples
- `helpfiles/<project>/bam_list.txt` rebuilt to include all BAM paths
- `helpfiles/<project>/hap_params.R` updated: add new sample names to `names_in_bam`
- `helpfiles/<project>/design.txt` updated: add rows for new samples
- New scan name so you don't overwrite the 3-rep results

**Step 1 — Align the new samples only:**

```bash
NEW_BARCODES=helpfiles/<project>/<project>_batch2.barcodes.txt
NN=$(wc -l < ${NEW_BARCODES})
sbatch --array=1-${NN} pipeline/scripts/fq2bam.sh \
    ${NEW_BARCODES} data/raw/<project>_batch2 data/bam/<project>
```

**Step 2 — Update your config files:**

- Rebuild `helpfiles/<project>/bam_list.txt` to include all BAMs (old + new)
- Update `hap_params.R`: add new sample names to `names_in_bam`
- Update `design.txt`: add rows for new samples

```bash
git add helpfiles/<project>/
git commit -m "add batch 2 samples to <project>"
git push
```

**Step 3 — Rerun from REFALT with a new scan name** (once alignment finishes):

```bash
bash pipeline/scripts/run_full_pipeline.sh \
    --skip-fq2bam \
    --project      <project> \
    --parfile      helpfiles/<project>/hap_params.R \
    --design       helpfiles/<project>/design.txt \
    --scan         <project>_6rep_smooth250 \
    --founders     A \
    --snp-table    pipeline/helpfiles/FREQ_SNPs_Apop.cM.txt.gz \
    --founder-list A1,A2,A3,A4,A5,A6,A7,AB8
```

Steps 3–4 must rerun with **all** BAMs because SNP calling and haplotype
inference are joint across all samples. Use a new scan name to preserve the
original results.

---

## Directory structure

```
XQTL2/                          ← this repo (pipeline — clone once per machine)
├── scripts/                    # Core pipeline scripts (tracked)
├── helpfiles/                  # Shared reference data (tracked)
│   ├── A_founders.bams.txt        # relative paths — no editing needed
│   ├── B_founders.bams.txt
│   ├── flymap.r6.txt
│   ├── het_bounds.txt
│   ├── FREQ_SNPs_Apop.cM.txt.gz
│   ├── FREQ_SNPs_Bpop.cM.txt.gz
│   ├── snp_tables/
│   └── generic_haplotype_parameters.R
├── ref/                        # Reference genome (not tracked — set up once)
└── data/founders/              # Founder BAMs (not tracked — set up once)

LongLab-XQTL/                   ← your project repo (any name)
├── pipeline -> ../XQTL2        # symlink — create after cloning
├── helpfiles/
│   └── <project>/              # your project config files (track in git)
│       ├── <project>.barcodes.txt
│       ├── bam_list.txt
│       ├── hap_params.R
│       └── design.txt
├── scripts_oneoffs/
│   └── <project>/
│       └── <project>_pipeline.sh
├── data/
│   ├── raw/<project>/          (not tracked)
│   └── bam/<project>/          (not tracked)
└── process/<project>/          (not tracked)
    ├── RefAlt.<chr>.txt
    ├── R.haps.<chr>.rds
    └── <scan_name>/
        ├── <scan_name>.scan.txt
        └── <scan_name>.tar.gz
```

To check what exists for a project on the cluster:

```bash
bash pipeline/scripts/show_project_layout.sh <project>
```
