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
6. Generate figures (`plot_pseudoscan.R`, `plot_H2_overlay.R`, `plot_freqsmooth_snp.R`)
7. Download results

**What's already in this repo:**

- All pipeline scripts (`scripts/`)
- Founder BAM paths (`helpfiles/founder.bams.txt`) — paths relative to this repo
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

The file `helpfiles/founder.bams.txt` uses `pipeline/data/founders/...` paths,
so it works directly from your project repo without any path manipulation.

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
to Drosophila; it can be adapted to any organism with pooled-sequencing XQTL data
by substituting the appropriate reference genome, genetic map, and founder BAMs.

### 4. Create your project repo

Your experimental data and project-specific scripts live in a separate repo
alongside XQTL2. Clone the [mylab-XQTL template](https://github.com/tdlong/mylab-XQTL)
and rename it to whatever suits your lab:

```bash
cd ..   # move up next to XQTL2/
git clone https://github.com/tdlong/mylab-XQTL.git LongLab-XQTL
cd LongLab-XQTL
ln -s ../XQTL2 pipeline   # all pipeline calls go through this symlink
```

Replace `LongLab-XQTL` with any name — `SmithLab-XQTL`, `MyProject`, whatever
makes sense for your group. The `pipeline` symlink is the only path that matters.

All your submission scripts call `pipeline/scripts/run_scan.sh` etc. When
XQTL2 is updated, run `git pull` inside the `/path/to/XQTL2` directory — your
scripts automatically use the new version via the symlink.

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

Save this file to `helpfiles/<project>/<project>.barcodes.txt`.

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

Create a file listing all BAM paths for your experiment — pooled samples plus
founders. Founders are pre-aligned and live in `pipeline/data/founders/`; the
provided `founder.bams.txt` and `B_founders.bams.txt` files already contain the
correct paths for your project repo.

For experiments using the standard A- or B-population founders, use the matching
file directly. AB8 is shared between populations. If your design crossed the
synthetic population to a tester strain or other reference genotype, treat that
strain as an additional "founder" and append its BAM path.

```bash
mkdir -p process/<project>
find data/bam/<project> -name "*.bam" -size +1G > helpfiles/<project>/bam_list.txt
cat pipeline/helpfiles/founder.bams.txt >> helpfiles/<project>/bam_list.txt
# for B-pop: cat pipeline/helpfiles/B_founders.bams.txt >> helpfiles/<project>/bam_list.txt

sbatch pipeline/scripts/bam2bcf2REFALT.sh \
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
Rscript pipeline/scripts/plot_pseudoscan.R \
    --scan   process/<project>/<scan_name>/<scan_name>.scan.txt \
    --out    process/<project>/<scan_name>/<scan_name>.wald.png \
    --format powerpoint \
    --threshold 10
```

Overlay two scans (e.g. male vs female):

```bash
Rscript pipeline/scripts/plot_pseudoscan.R \
    --scan   process/<project>/<scan_M>/<scan_M>.scan.txt \
    --scan   process/<project>/<scan_F>/<scan_F>.scan.txt \
    --label  Male --label Female \
    --colour "#1F78B4" --colour "#E31A1C" \
    --out    process/<project>/MF_overlay.png \
    --format powerpoint --threshold 10
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

`plot_pseudoscan.R` and `plot_freqsmooth_snp.R` also accept `--label` and
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
bash scripts_oneoffs/malathion/malathion_pipeline.sh
```

The script chains Steps 3–6 with SLURM dependencies and prints each job ID.
When complete, download the results tarball and check for a signal on chr3R
around 9 Mb (the *Ace* locus — a known malathion resistance gene).

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

```bash
#!/bin/bash
set -e

PROJECT=myproject
PARFILE=helpfiles/${PROJECT}/hap_params.R
DESIGN=helpfiles/${PROJECT}/design.txt
SCAN=${PROJECT}_6rep_smooth250
SNP_TABLE=pipeline/helpfiles/FREQ_SNPs_Apop.cM.txt.gz
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8

# Align NEW samples only
NEW_BARCODES=helpfiles/${PROJECT}/${PROJECT}_batch2.barcodes.txt
NN=$(wc -l < ${NEW_BARCODES})
jid_bam=$(sbatch --parsable --array=1-${NN} pipeline/scripts/fq2bam.sh \
    ${NEW_BARCODES} data/raw/${PROJECT}_batch2 data/bam/${PROJECT})

# Rebuild bam_list with old + new, rerun REFALT
find data/bam/${PROJECT} -name "*.bam" -size +1G > helpfiles/${PROJECT}/bam_list.txt
grep "A" pipeline/helpfiles/founder.bams.txt | sed 's|^|pipeline/|' \
    >> helpfiles/${PROJECT}/bam_list.txt

jid_refalt=$(sbatch --parsable --dependency=afterok:${jid_bam} \
    pipeline/scripts/bam2bcf2REFALT.sh \
    helpfiles/${PROJECT}/bam_list.txt process/${PROJECT})

# Rerun haplotypes (hap_params.R must be updated with all sample names)
jid_haps=$(sbatch --parsable --dependency=afterok:${jid_refalt} \
    --array=1-5 pipeline/scripts/REFALT2haps.sh \
    --parfile ${PARFILE} --dir process/${PROJECT})

# Scan
scan_out=$(bash pipeline/scripts/run_scan.sh \
    --design ${DESIGN} --dir process/${PROJECT} \
    --scan ${SCAN} --after ${jid_haps})
jid_hap=$(echo "$scan_out" | grep "^done:" | awk '{print $2}')

snp_out=$(bash pipeline/scripts/run_snp_scan.sh \
    --design ${DESIGN} --dir process/${PROJECT} \
    --scan ${SCAN} --snp-table ${SNP_TABLE} --founders ${FOUNDERS})
jid_snp=$(echo "$snp_out" | grep "^done:" | awk '{print $2}')

# Figures
SCAN_DIR=process/${PROJECT}/${SCAN}
sbatch --dependency=afterok:${jid_hap},afterok:${jid_snp} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
    --wrap="module load R/4.2.2 && \
Rscript pipeline/scripts/plot_pseudoscan.R \
    --scan ${SCAN_DIR}/${SCAN}.scan.txt \
    --out  ${SCAN_DIR}/${SCAN}.wald.png --format powerpoint --threshold 10 && \
Rscript pipeline/scripts/plot_H2_overlay.R \
    --scan ${SCAN_DIR}/${SCAN}.scan.txt \
    --out  ${SCAN_DIR}/${SCAN}.H2.png --format powerpoint && \
Rscript pipeline/scripts/plot_freqsmooth_snp.R \
    --scan ${SCAN_DIR}/${SCAN}.snp_scan.txt \
    --out  ${SCAN_DIR}/${SCAN}.snp.wald.png --format powerpoint --threshold 10 && \
cd ${SCAN_DIR} && tar -czf ${SCAN}.tar.gz *.txt *.png"
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
│   ├── founder.bams.txt        # relative paths — no editing needed
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
