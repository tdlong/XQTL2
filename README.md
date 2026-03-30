# XQTL Pipeline

## Overview

This pipeline takes pooled-sequencing XQTL data from raw reads to genome scans
and publication figures. Everything runs on a SLURM cluster from the project
root directory. You submit one script and walk away.

**Pipeline at a glance:**

1. Get raw reads
2. Align reads (`fq2bam.sh`)
3. Generate allele counts (`bam2bcf2REFALT.sh`)
4. Call haplotypes (`REFALT2haps.sh`)
5. Scan
   - 5a. Haplotype scan (`run_scan.sh` — smooth, Wald test, H², concat)
   - 5b. SNP scan (`run_snp_scan.sh` — optional, imputed SNP-level test)
6. Generate figures (Rscript commands, submitted via sbatch)
7. Download results

**What's already in this repo:**

- All pipeline scripts (`scripts/`)
- Founder bam paths (`helpfiles/founder.bams.txt`) — pre-aligned, shared on cluster
- SNP frequency tables (`helpfiles/FREQ_SNPs_Apop.cM.txt.gz`, `FREQ_SNPs_Bpop.cM.txt.gz`)
- Genetic map (`helpfiles/flymap.r6.txt`)
- A complete worked example (`helpfiles/malathion_test/`)

**What you need to provide per experiment:**

- Raw sequencing reads (from the sequencing core)
- A barcode file mapping barcodes → sample names (Step 2)
- A haplotype parameters file listing your founders and samples (Step 4)
- A design file describing your experimental layout (Step 5)

**What you get at the end:**

- `*.scan.txt` — genome-wide haplotype scan (Wald -log10p, Falconer H², Cutler H²)
- `*.snp_scan.txt` — SNP-level Wald test at every imputed SNP (optional)
- `*.meansBySample.txt` — smoothed founder frequencies per sample (QC)
- Manhattan plots and heritability figures (PNG)
- Tarballs of everything, ready to scp down

**Note:** SLURM headers use `-A tdlong_lab -p standard`. If you're in a
different lab, change the account/partition in the SBATCH scripts.

A legacy scan without smoothing (`haps2scan.Apr2025.sh`) is also available.

### SLURM resource requirements

The cluster's standard partition provides max 6 GB per core; highmem provides
10 GB per core. Always use `--mem-per-cpu` (not `--mem`). See `Slurm.md` for
full partition details.

Every script explicitly requests memory. Scan steps were profiled with `seff`
on the malathion test dataset (2 replicates, 4 samples). Larger experiments
will scale proportionally. Baseline request is 1 CPU × 3G — there is little
gain in requesting less. For Steps 5–6, `run_scan.sh` and `run_snp_scan.sh`
accept `--mem-per-cpu`, `--cpus-per-task`, `-p`, and `-A` to override
resources and partition for all jobs they submit.

| Script | Step | Partition | CPUs | Mem/CPU | Time | Profiled (malathion test) |
|--------|------|-----------|------|---------|------|--------------------------|
| `fq2bam.sh` | 2 | standard | 4 | 6G | 1 day | `bwa -t 4` uses 4 threads; `java -Xmx20g` needs ~20G |
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
Only include founders for your population — grep `"A"` for A-pop or `"B"` for B-pop
(AB8 is shared and matched by both).

```bash
mkdir -p process/<project>
find data/bam/<project> -name "*.bam" -size +1G > helpfiles/<project>/bam_list.txt
grep "A" helpfiles/founder.bams.txt >> helpfiles/<project>/bam_list.txt   # or "B" for B-pop

sbatch scripts/bam2bcf2REFALT.sh \
    helpfiles/<project>/bam_list.txt \
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
Wald test and heritability estimates, and concatenating chromosomes. One command
per scan.

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
| `--mem-per-cpu` | 3G | Memory per CPU for all jobs (SLURM `--mem-per-cpu`) |
| `--cpus-per-task` | 1 | CPUs for all jobs (SLURM `--cpus-per-task`) |
| `-p` / `--partition` | standard | SLURM partition (e.g. `highmem`) |
| `-A` / `--account` | tdlong_lab | SLURM account to charge |
| `--after` | (none) | SLURM job ID to wait on before starting |

### Smoothing window

The default smoothing window is 250 kb. This was chosen by comparing 125 kb and
250 kb on real data (malathion experiment): plotting smoothed founder haplotype
frequencies versus genomic position and comparing against expectations from
simulations. At 250 kb the founder frequency estimates are stable without
over-smoothing genuine biological signal. You can override with `--smooth 125`
or any other value.

### Unresolvable founder masking

When two or more founders have nearly identical haplotypes in a region, the
haplotype estimator cannot resolve their individual frequencies — only their sum
is constrained. The smoothing step detects these cases using the `Groups` output
from REFALT2haps and masks the unresolvable founders to NA before smoothing.
The NA-safe running mean then interpolates from flanking valid data.

This works well for small gaps where only some founders are ambiguous. For large
regions where **all** founders are simultaneously unresolvable (e.g., chrX
20.5–22.6 Mb in some datasets), the smoothing kernel cannot reach valid anchors
and those windows will be absent from the scan output. Any QTL in such a region
is undetectable regardless of method — there is no founder-level signal to recover.

### What run_scan.sh submits

For reference, `run_scan.sh` chains these SLURM jobs automatically:

1. `smooth_haps.sh` — smooth haplotype frequencies and covariances (5 array tasks)
2. `hap_scan.sh` — Wald test + heritability at each haplotype window (5 array tasks, after #1)
3. `concat_scans.sh` — merge per-chromosome files and generate Manhattan plots (after #2)

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

The SNP scan imputes per-SNP ALT allele frequencies from the smoothed haplotype
estimates produced in Step 5a, then runs a Wald test (df=1) at every SNP.
This tests at individual SNP positions rather than haplotype windows, but
the signal comes from the same smoothed haplotype estimates — it is not
independent of the haplotype scan.

Outputs: `snp_scan.txt` (Wald -log10(p) per SNP) and `snp_meansBySample.txt`
(imputed ALT frequency per SNP per treatment per replicate). Heritability is
not estimated at the SNP level — H² is a property of genomic regions, not
individual SNPs (see Step 5a for H² estimates).

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
    --founders  A1,A2,A3,A4,A5,A6,A7,AB8
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
| `--mem-per-cpu` | 3G | Memory per CPU for all jobs (SLURM `--mem-per-cpu`) |
| `--cpus-per-task` | 1 | CPUs for all jobs (SLURM `--cpus-per-task`) |
| `-p` / `--partition` | standard | SLURM partition (e.g. `highmem`) |
| `-A` / `--account` | tdlong_lab | SLURM account to charge |
| `--after` | (none) | SLURM job ID to wait on (e.g. from run_scan.sh) |

---

## Step 6 — Generate publication figures

Three plotting scripts produce 5-panel (per-chromosome) figures. Run these
on the cluster after the scan completes — either submit via sbatch (see worked
example) or run interactively. Each takes the scan output file as a
command-line argument.

### Haplotype Wald scan

```bash
Rscript scripts/plot_pseudoscan.R \
    --scan   process/<project>/<scan_name>/<scan_name>.scan.txt \
    --out    process/<project>/<scan_name>/wald.png \
    --format powerpoint \
    --threshold 10
```

Overlay two scans (e.g. male vs female):

```bash
Rscript scripts/plot_pseudoscan.R \
    --scan   process/<project>/<scan_M>/<scan_M>.scan.txt \
    --scan   process/<project>/<scan_F>/<scan_F>.scan.txt \
    --label  Male --label Female \
    --colour "#1F78B4" --colour "#E31A1C" \
    --out    process/<project>/MF_overlay.png \
    --format powerpoint --threshold 10
```

### Heritability overlay (Falconer + Cutler)

```bash
Rscript scripts/plot_H2_overlay.R \
    --scan   process/<project>/<scan_name>/<scan_name>.scan.txt \
    --out    process/<project>/<scan_name>/H2.png \
    --format powerpoint
```

### SNP scan

```bash
Rscript scripts/plot_freqsmooth_snp.R \
    --scan   process/<project>/<scan_name>/<scan_name>.snp_scan.txt \
    --out    process/<project>/<scan_name>/snp_wald.png \
    --format powerpoint --threshold 10
```

### Common options

All three scripts accept these arguments:

| Flag | Description |
|------|-------------|
| `--scan <file>` | Input scan file (required; repeat for overlays in pseudoscan/snp) |
| `--out <file>` | Output PNG path (required) |
| `--format <name>` | Size/DPI preset (default: `powerpoint`) |
| `--threshold <n>` | Dashed horizontal line at this y value |
| `--genes <file>` | Tab-delimited gene annotations (columns: `name`, `chr`, `pos_mb`) |
| `--peaks <file>` | Tab-delimited peak annotations (columns: `label`, `chr`, `pos_mb`) |
| `--height <in>` | Override figure height in inches (default: 1.4 per chromosome) |

`plot_pseudoscan.R` and `plot_freqsmooth_snp.R` also accept `--label` and
`--colour` (one per `--scan`, for overlays).

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

### Gene and peak annotation files

To label genes or peaks on any figure, create a tab-delimited text file and
pass it with `--genes` or `--peaks`. These work with all three plot engines.

`helpfiles/<project>/genes.txt`:

```
name	chr	pos_mb
Ace	chr3R	9.07
Cyp6g1	chr2R	12.19
```

`helpfiles/<project>/peaks.txt`:

```
label	chr	pos_mb
peak1	chr3R	9.1
```

```bash
Rscript scripts/plot_pseudoscan.R \
    --scan ... --out ... --format powerpoint \
    --genes helpfiles/<project>/genes.txt \
    --peaks helpfiles/<project>/peaks.txt
```

---

## Step 7 — Download results

The figure step (Step 6) bundles everything into a single tarball at the end.

```bash
scp <user>@<cluster>:<project_path>/process/<project>/<scan_name>/<scan_name>.tar.gz .
tar xzf <scan_name>.tar.gz
```

**What's in the tarball** (all `.txt` and `.png` files in the scan directory):

| File | Contents |
|------|----------|
| `<scan>.scan.txt` | Haplotype scan results (one row per window) |
| `<scan>.meansBySample.txt` | Smoothed founder haplotype frequencies (one row per window × treatment × replicate × founder) |
| `<scan>.snp_scan.txt` | SNP scan results (one row per SNP; if SNP scan was run) |
| `<scan>.snp_meansBySample.txt` | Imputed SNP ALT frequencies (one row per SNP × treatment × replicate; if SNP scan was run) |
| `<scan>.wald.png` | 5-panel haplotype Wald Manhattan |
| `<scan>.H2.png` | 5-panel Falconer + Cutler heritability overlay |
| `<scan>.snp.wald.png` | 5-panel SNP Wald Manhattan (if SNP scan was run) |
| `<scan>.5panel.*.png`, `<scan>.Manhattan.png` | Quick-look plots from concat step |

### Output file formats

**`<scan>.scan.txt`** — haplotype scan (produced by `hap_scan.R`, one row per haplotype window):

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `pos` | Window center (bp) |
| `Wald_log10p` | -log10(p) from Wald test |
| `Falc_H2` | Falconer heritability estimate |
| `Cutl_H2` | Cutler heritability estimate |
| `cM` | Genetic map position (centiMorgans) |

**`<scan>.meansBySample.txt`** — smoothed founder frequencies (produced by `smooth_haps.R`):

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `pos` | Window center (bp) |
| `TRT` | Treatment: `C` (control) or `Z` (selected) |
| `REP` | Replicate number |
| `founder` | Founder name (e.g. A1, A2, ..., AB8) |
| `freq` | Smoothed founder haplotype frequency |

**`<scan>.snp_scan.txt`** — SNP scan (produced by `snp_scan.R`, one row per SNP):

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `pos` | SNP position (bp) |
| `Wald_log10p` | -log10(p) from Wald test |
| `cM` | Genetic map position (centiMorgans) |
| `n_informative_founders` | Number of founders carrying the ALT allele |

**`<scan>.snp_meansBySample.txt`** — imputed SNP frequencies (produced by `snp_scan.R`):

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `pos` | SNP position (bp) |
| `TRT` | Treatment: `C` (control) or `Z` (selected) |
| `REP` | Replicate number |
| `F_alt` | Imputed ALT allele frequency |
| `cM` | Genetic map position (centiMorgans) |

### Interactive exploration (optional)

`scripts/XQTL_plotting_functions.R` provides functions for zooming into peaks
and making regional plots in an R session:

```r
source("scripts/XQTL_plotting_functions.R")
df1 <- as_tibble(read.table("SCAN_NAME.scan.txt"))
df2 <- as_tibble(read.table("SCAN_NAME.meansBySample.txt"))

XQTL_Manhattan_5panel(df1, cM = FALSE)
XQTL_region(df1, "chr3R", 18250000, 19000000, "Wald_log10p")
XQTL_change_average(df2, "chr3R", 18250000, 19000000)
XQTL_combined_plot(df1, df2, "chr3R", 18250000, 19000000)
```

---

## Worked example — end-to-end pipeline

Copy this to `scripts_oneoffs/<project>_pipeline.sh`, fill in the variables at
the top, and run it. Everything runs on the cluster — scans, figures, tarballs.
When it's done, scp the results down.

**What you need before running:** the four input files from Steps 2–5 above
(barcode file, hap_params.R, design file) and raw reads in `data/raw/<project>/`.

**What you get when it finishes:** for each scan, a directory containing
haplotype scan tables, SNP scan tables, Manhattan plots, heritability figures,
and tarballs ready to download.

```bash
#!/bin/bash
set -e

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  Fill these in for your experiment                                       ║
# ╚═══════════════════════════════════════════════════════════════════════════╝
PROJECT=myproject
BARCODES=helpfiles/${PROJECT}/${PROJECT}.barcodes.txt
PARFILE=helpfiles/${PROJECT}/hap_params.R
DESIGN=helpfiles/${PROJECT}/design.txt
SCAN=${PROJECT}_smooth250
SNP_TABLE=helpfiles/FREQ_SNPs_Apop.cM.txt.gz   # or Bpop for B-population
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8

# ── Step 2: Align reads ─────────────────────────────────────────────────────
NN=$(wc -l < ${BARCODES})
mkdir -p data/bam/${PROJECT}
jid_bam=$(sbatch --parsable --array=1-${NN} scripts/fq2bam.sh \
    ${BARCODES} data/raw/${PROJECT} data/bam/${PROJECT})

# ── Step 3: REFALT counts ───────────────────────────────────────────────────
mkdir -p process/${PROJECT}
find data/bam/${PROJECT} -name "*.bam" -size +1G > helpfiles/${PROJECT}/bam_list.txt
grep "A" helpfiles/founder.bams.txt >> helpfiles/${PROJECT}/bam_list.txt   # or "B" for B-pop
jid_refalt=$(sbatch --parsable --dependency=afterok:${jid_bam} \
    scripts/bam2bcf2REFALT.sh helpfiles/${PROJECT}/bam_list.txt process/${PROJECT})

# ── Step 4: Call haplotypes ──────────────────────────────────────────────────
jid_haps=$(sbatch --parsable --dependency=afterok:${jid_refalt} \
    --array=1-5 scripts/REFALT2haps.sh \
    --parfile ${PARFILE} --dir process/${PROJECT})

# ── Step 5a: Haplotype scan (smooth → Wald test + H² → concat) ──────────────
scan_out=$(bash scripts/run_scan.sh \
    --design ${DESIGN} \
    --dir    process/${PROJECT} \
    --scan   ${SCAN} \
    --after  ${jid_haps})
echo "$scan_out"
jid_hap=$(echo "$scan_out" | grep "^done:" | awk '{print $2}')

# ── Step 5b: SNP scan (optional — delete this block if not needed) ──────────
snp_out=$(bash scripts/run_snp_scan.sh \
    --design    ${DESIGN} \
    --dir       process/${PROJECT} \
    --scan      ${SCAN} \
    --snp-table ${SNP_TABLE} \
    --founders  ${FOUNDERS})
echo "$snp_out"
jid_snp=$(echo "$snp_out" | grep "^done:" | awk '{print $2}')

# ── Step 6: Figures + final tarball (runs after all scans finish) ─────────────
SCAN_DIR=process/${PROJECT}/${SCAN}
sbatch --dependency=afterok:${jid_hap},afterok:${jid_snp} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
    --wrap="module load R/4.2.2 && \
Rscript scripts/plot_pseudoscan.R \
    --scan      ${SCAN_DIR}/${SCAN}.scan.txt \
    --out       ${SCAN_DIR}/${SCAN}.wald.png \
    --format    powerpoint \
    --threshold 10 && \
Rscript scripts/plot_H2_overlay.R \
    --scan   ${SCAN_DIR}/${SCAN}.scan.txt \
    --out    ${SCAN_DIR}/${SCAN}.H2.png \
    --format powerpoint && \
Rscript scripts/plot_freqsmooth_snp.R \
    --scan      ${SCAN_DIR}/${SCAN}.snp_scan.txt \
    --out       ${SCAN_DIR}/${SCAN}.snp.wald.png \
    --format    powerpoint \
    --threshold 10 && \
cd ${SCAN_DIR} && tar -czf ${SCAN}.tar.gz *.txt *.png"

echo "All jobs submitted."
echo "When done: scp <user>@<cluster>:$(pwd)/${SCAN_DIR}/${SCAN}.tar.gz ."
```

---

## Worked example — adding replicates to an existing experiment

You sequenced 3 replicates, ran the pipeline, then sequenced 3 more. Now you
want to reanalyze with all 6. The existing bams and process directory are still
on the cluster — you only need to align the new samples, then rerun from Step 3
onward with all bams combined.

**What changes:**
- New barcode file for the new samples only (or append to original)
- `helpfiles/<project>/bam_list.txt` rebuilt to include old + new bam paths
- `helpfiles/<project>/hap_params.R` updated: add new sample names to `names_in_bam`
- `helpfiles/<project>/design.txt` updated: add rows for new samples
- New scan name (e.g. `myproject_6rep_smooth250`) so you don't overwrite the 3-rep results

**What stays the same:** founder bams, reference genome, SNP table, all scripts.

```bash
#!/bin/bash
set -e

PROJECT=myproject
PARFILE=helpfiles/${PROJECT}/hap_params.R
DESIGN=helpfiles/${PROJECT}/design.txt
SCAN=${PROJECT}_6rep_smooth250
SNP_TABLE=helpfiles/FREQ_SNPs_Apop.cM.txt.gz
FOUNDERS=A1,A2,A3,A4,A5,A6,A7,AB8

# ── Step 2: Align NEW samples only ───────────────────────────────────────────
NEW_BARCODES=helpfiles/${PROJECT}/${PROJECT}_batch2.barcodes.txt
NN=$(wc -l < ${NEW_BARCODES})
jid_bam=$(sbatch --parsable --array=1-${NN} scripts/fq2bam.sh \
    ${NEW_BARCODES} data/raw/${PROJECT}_batch2 data/bam/${PROJECT})

# ── Step 3: Rebuild bams list (old + new) and rerun REFALT ───────────────────
#   Old bams are already in data/bam/<project>/ from the first run.
#   After alignment finishes, combine all bam paths into one file.
find data/bam/${PROJECT} -name "*.bam" -size +1G > helpfiles/${PROJECT}/bam_list.txt
grep "A" helpfiles/founder.bams.txt >> helpfiles/${PROJECT}/bam_list.txt   # or "B" for B-pop

jid_refalt=$(sbatch --parsable --dependency=afterok:${jid_bam} \
    scripts/bam2bcf2REFALT.sh helpfiles/${PROJECT}/bam_list.txt process/${PROJECT})

# ── Step 4: Rerun haplotype calling with all samples ─────────────────────────
#   Make sure hap_params.R has been updated with all sample names.
jid_haps=$(sbatch --parsable --dependency=afterok:${jid_refalt} \
    --array=1-5 scripts/REFALT2haps.sh \
    --parfile ${PARFILE} --dir process/${PROJECT})

# ── Step 5a: Haplotype scan ─────────────────────────────────────────────────
scan_out=$(bash scripts/run_scan.sh \
    --design ${DESIGN} \
    --dir    process/${PROJECT} \
    --scan   ${SCAN} \
    --after  ${jid_haps})
echo "$scan_out"
jid_hap=$(echo "$scan_out" | grep "^done:" | awk '{print $2}')

# ── Step 5b: SNP scan ───────────────────────────────────────────────────────
snp_out=$(bash scripts/run_snp_scan.sh \
    --design    ${DESIGN} \
    --dir       process/${PROJECT} \
    --scan      ${SCAN} \
    --snp-table ${SNP_TABLE} \
    --founders  ${FOUNDERS})
echo "$snp_out"
jid_snp=$(echo "$snp_out" | grep "^done:" | awk '{print $2}')

# ── Step 6: Figures + tarball ────────────────────────────────────────────────
SCAN_DIR=process/${PROJECT}/${SCAN}
sbatch --dependency=afterok:${jid_hap},afterok:${jid_snp} \
    -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=3G --time=1:00:00 \
    --wrap="module load R/4.2.2 && \
Rscript scripts/plot_pseudoscan.R \
    --scan      ${SCAN_DIR}/${SCAN}.scan.txt \
    --out       ${SCAN_DIR}/${SCAN}.wald.png \
    --format    powerpoint \
    --threshold 10 && \
Rscript scripts/plot_H2_overlay.R \
    --scan   ${SCAN_DIR}/${SCAN}.scan.txt \
    --out    ${SCAN_DIR}/${SCAN}.H2.png \
    --format powerpoint && \
Rscript scripts/plot_freqsmooth_snp.R \
    --scan      ${SCAN_DIR}/${SCAN}.snp_scan.txt \
    --out       ${SCAN_DIR}/${SCAN}.snp.wald.png \
    --format    powerpoint \
    --threshold 10 && \
cd ${SCAN_DIR} && tar -czf ${SCAN}.tar.gz *.txt *.png"

echo "All jobs submitted."
```

The key difference from a fresh run: you only align the new samples (Step 2),
but Steps 3–4 must rerun with **all** bams because SNP calling and haplotype
inference are joint across all samples. Use a new scan name so the original
results are preserved for comparison.

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
│       ├── bam_list.txt                  (Step 3)
│       ├── hap_params.R                  (Step 4)
│       └── design.txt                    (Step 5)
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
