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
3. Generate allele counts (`build_catalog.sh` once, then `call_samples.sh`)
   - Per-sample QC (`refalt_qc.R`, `bam_qc.sh`) — automatic in `run_full_pipeline.sh`
4. Call haplotypes (`REFALT2haps.sh`)
5. Scan
   - 5a. Haplotype scan (`run_scan.sh` — smooth, Wald test, H², concat)
   - 5b. SNP scan (`run_snp_scan.sh` — optional, imputed SNP-level test)
6. Generate figures (`plot_5panel.R`, `plot_manhattan.R`, `plot_H2_overlay.R`, `plot_freqsmooth_snp.R`)
7. Download results

**What's already in this repo:**

- All pipeline scripts (`scripts/`)
- Founder BAM paths (`helpfiles/A_founders.bams.txt`) — paths relative to this repo
- Physical-to-genetic map (`helpfiles/flymap.r6.txt`)
- Heterochromatic boundary definitions (`helpfiles/het_bounds.txt`)
- Generic haplotype parameter template (`helpfiles/generic_haplotype_parameters.R`)

**What you need to provide per experiment:**

- Raw sequencing reads (from the sequencing core)
- A barcode file mapping barcodes → sample names (Step 2)
- A haplotype parameters file listing your founders, samples, and window sizes (Step 4)
- A design file describing your experimental layout (Step 5)

**What you get at the end:**

- `*.scan.txt` — genome-wide haplotype scan (Wald -log10p, H²)
- `*.snp_scan.txt` — SNP-level Wald test at every imputed SNP (optional)
- `Calls/refalt_qc.txt`, `Calls/bam_qc.txt` — per-sample coverage and alignment QC
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
| `catalog_build.sh` (×5 chr) | 3 | standard | 1 | 6G | ~4 hr | founder calling; once per catalog |
| `catalog_count.sh` (×N samples) | 3 | standard | 1 | 6G | 2–11 min | one job per sample at a fixed catalog |
| `REFALT2haps.sh` | 4 | highmem | 1 | 10G | 1 day | large haplotype matrices require highmem |
| `smooth_haps.sh` | 5a | highmem | 1 | 10G | 4 hr | covariance uncount(64) needs ~8G on large chromosomes |
| `smooth_r2_diag.sh` | 5a | standard | 1 | 4G | 30 min | R² calibration; runs between smooth and hap_scan |
| `hap_scan.sh` | 5a | standard | 1 | 3G | 4 hr | 307 MB / 5:12 wall; reads R² from smooth_r2.txt |
| `refalt_qc.R` + `bam_qc.sh` | 3 | standard | 1 | 8G | 2 hr | one pass over RefAlt + one flagstat per BAM |
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

The default SNP caller is the **founder-catalog caller**: decide the SNP list
once from the founders, then count each sample's REF/ALT reads at those fixed
positions. It is two deliberate commands — `build_catalog.sh` (3c) then
`call_samples.sh` (3d) — and it produces `Calls/RefAlt.<chr>.txt`, the sole input
to everything downstream.

Why it works this way, the exact filters, resource sizing, and what each kind of
rerun costs are in the **"Founder-catalog caller"** appendix at the end of this
README. The superseded joint QUAL caller is in
[`scripts/legacy/`](scripts/legacy/README.md), kept only for reproducing older
analyses.

### 3a — Choose your founders

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

### 3b — Build the bam list

Create `helpfiles/<project>/bam_list.txt` — one BAM path per line: your sample
BAMs, then whichever founders you want as columns in `RefAlt`. `call_samples.sh`
counts **exactly** the BAMs you list, so this file is the record of what went into
the callset. Build a draft from your BAM directory, append founders, then
**review it** before submitting.

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

### 3c — Build the catalog (once per population)

```bash
mkdir -p process/<project>
bash pipeline/scripts/build_catalog.sh \
    --founders pipeline/helpfiles/A_founders.bams.txt \
    --out      process/<project>/Catalog
```

Founder calling runs as a 5-task array (one chromosome each: chrX, chr2L, chr2R,
chr3L, chr3R), then gathers and applies thresholds. This is the slow part (~4 h
per chromosome) and you do it **once** — every later run, and every added sample,
reuses it.

Building **overwrites** `--out` and has no staleness check, so it is a deliberate
act: nothing else in the pipeline builds a catalog for you. Thresholds
(`--min-dp`, `--maxaf`, `--snpgap`, `--exempt-founders`) can be re-cut in seconds
afterwards without recalling the founders — see the appendix.

### 3d — Count your samples against it

```bash
JID_REFALT=$(bash pipeline/scripts/call_samples.sh \
    --catalog process/<project>/Catalog \
    --bamlist helpfiles/<project>/bam_list.txt \
    --dir     process/<project>)
```

One job per sample (2–11 min each), then a merge into
`process/<project>/Calls/RefAlt.<chr>.txt`. It prints the merge job ID, so you can
chain Step 4 on it with `--after ${JID_REFALT}`.

**Adding samples later costs only the new samples.** Rerun `call_samples.sh` with a
`--bamlist` of *just the new BAMs*: the catalog is reused, existing per-sample
counts are left alone, and `RefAlt.<chr>.txt` is re-merged to include everything.
See **"Worked example — adding replicates to an existing experiment"**.

### Check the samples before calling haplotypes

A sample with poor coverage does not fail loudly. `est_hap2` fits each sample on
its own, so a bad BAM cannot corrupt its neighbours — but a thin one returns
haplotype frequencies that are badly wrong and still satisfy every constraint the
pipeline checks. The Wald test partly protects itself (a large reconstruction
covariance drives that replicate's effective N down); H² does not use the
reconstruction covariance at all, taking its variance from the replicates.

```bash
Rscript pipeline/scripts/refalt_qc.R \
    --dir     process/<project> \
    --parfile helpfiles/<project>/hap_params.R
```

**`run_full_pipeline.sh` runs this for you**, as soon as RefAlt exists — it is
submitted alongside the haplotype step rather than in the chain, so it reports but
never gates. Run it by hand as above when driving the steps individually, or to
re-check an existing project.

Writes `Calls/refalt_qc.txt`, one row per sample per chromosome, and prints
anything flagged. Every metric is that sample's own — nothing is relative to
other samples.

| column | reads on |
|---|---|
| `median_depth` | library success. Median, since collapsed repeats inflate the mean |
| `pct_zero`, `pct_lt5` | sites with no reads, and with too few to be informative |
| `sites_per_window` | covered catalog sites in a typical haplotype window |
| `pct_window_covered` | that as a % of what the catalog offers — **the one flagged on**, since it separates a thin sample from a sparse catalog |
| `disp_k` | overdispersion of depth relative to Poisson |
| `flag` | `OK` / `LOW_COVERAGE` / `PATCHY` |

`disp_k` catches the sample that has plenty of reads piled into too few places.
It is a negative-binomial size fitted by moments over depths within [0.25×, 4×]
of the median — the band matters, because a collapsed-repeat tail of 1% of sites
at 20× depth drags an untrimmed estimate from a true 10 down to 0.35. These
libraries are never Poisson; the number is for ranking samples against each
other, not an absolute. It is deliberately depth-invariant, unlike `var/mean`,
which grows with mean depth and so cannot compare samples sequenced to different
depths.

For alignment-side numbers — whether a thin sample had few reads or plenty that
did not map:

```bash
bash pipeline/scripts/bam_qc.sh --bamlist helpfiles/<project>/bam_list.txt \
    --out process/<project>/Calls/bam_qc.txt
```

Joins to `refalt_qc.txt` on `sample`. `has_unmapped` is reported because
`pct_mapped` only means something if unmapped reads are still present — a BAM
that arrived via `--skip-fq2bam` may have had them stripped, in which case the
percentage is a meaningless 100. There is no duplicate metric: these libraries
are not deduplicated by design, since at 50–200× reads sharing a start site occur
legitimately, especially from Nextera.

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

```bash
JID_HAPS=$(bash pipeline/scripts/run_haps.sh \
    --after   $JID_REFALT \
    --parfile helpfiles/<project>/hap_params.R \
    --dir     process/<project>)
```

This submits a 5-task array (one chromosome per task) and prints the job ID.
The `--after` flag makes it wait for REFALT to finish before starting.

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
| `--sex` | mixed | Sex of the pools in this scan: `M`, `F` or `mixed`. chrX only — see below |
| `--mem-per-cpu` | 3G | Memory per CPU for all jobs |
| `--cpus-per-task` | 1 | CPUs for all jobs |
| `-p` / `--partition` | standard | SLURM partition |
| `-A` / `--account` | tdlong_lab | SLURM account |
| `--after` | (none) | SLURM job ID to wait on (e.g. from REFALT2haps) |

#### `--sex` — required for single-sex experiments

`Num` counts flies, and the scan works in chromosomes: `2 * Num` of them. That
is right on an autosome, where every fly carries two. On chrX it depends on the
sex of the flies, so `Num` is scaled by half the number of X chromosomes per
fly:

| `--sex` | Scaling | Chromosomes |
|---------|---------|-------------|
| `mixed` (default) | 0.75 (1.5 X per fly) | `1.5 * Num` |
| `F` | 1.0 (2 X per fly) | `2.0 * Num` |
| `M` | 0.5 (1 X per fly) | `1.0 * Num` |

The default `mixed` is what the pipeline has always assumed, so existing
projects are unaffected. If your flies are single-sex you must set it: at
`mixed` a male pool is credited with 1.5× the X chromosomes it carries and a
female pool with 0.75×, which inflates chrX statistics for males and deflates
them for females. Autosomes are never affected.

The flag is per scan, not per sample, because a scan contrasting a male pool
against a female one would confound the treatment effect with sex. If your
experiment has both, run one scan per sex and compare the two scans:

```bash
bash pipeline/scripts/run_scan.sh \
    --design helpfiles/<project>/design.males.txt \
    --dir    process/<project> \
    --scan   <project>_males_smooth250 \
    --sex    M

bash pipeline/scripts/run_scan.sh \
    --design helpfiles/<project>/design.females.txt \
    --dir    process/<project> \
    --scan   <project>_females_smooth250 \
    --sex    F
```

Each design file lists only that sex's pools, and each scan gets its own
output directory, so the two never mix.

### Option B — Step by step

For reference or debugging, `run_scan.sh` internally chains these four jobs:

**1. Smooth haplotype frequencies** (5-task array, one chromosome per task):

```bash
sbatch --array=1-5 pipeline/scripts/smooth_haps.sh \
    --rfile     helpfiles/<project>/design.txt \
    --dir       process/<project> \
    --outdir    <scan_name> \
    --smooth-kb 250 \
    --sex       mixed
```

**2. R² smoothing diagnostic** — computes the R² between smoothed and raw
haplotype frequencies across all founders and pools, writes
`<scan_name>.smooth_r2.txt` to the scan directory. The haplotype scan
reads this file to apply the calibration correction automatically.

```bash
sbatch pipeline/scripts/smooth_r2_diag.sh \
    --hapsdir   process/<project>/Haps \
    --smoothdir process/<project>/Scans/<scan_name> \
    --scan      <scan_name> \
    --rfile     helpfiles/<project>/design.txt
```

**3. Haplotype scan** — Wald test + heritability at each window (5-task array,
after smooth_r2_diag):

```bash
sbatch --array=1-5 pipeline/scripts/hap_scan.sh \
    --rfile  helpfiles/<project>/design.txt \
    --dir    process/<project> \
    --outdir <scan_name>
```

**4. Concatenate chromosomes and generate quick-look plots:**

```bash
bash pipeline/scripts/concat_scans.sh process/<project>/Scans/<scan_name>
```

### Smoothing window

The default smoothing window is 250 kb, chosen by comparing 125 kb and 250 kb
on real data (malathion experiment). At 250 kb the founder frequency estimates
are stable without over-smoothing genuine biological signal. Override with
`--smooth 125` or any value.

### Legacy scan (no smoothing)

An older scan method without haplotype smoothing is preserved in
[`scripts/legacy/`](scripts/legacy/README.md) for reproducing pre-smoothing
results. It is not part of the current pipeline:

```bash
sbatch --array=1-5 pipeline/scripts/legacy/haps2scan.Apr2025.sh \
    --rfile  helpfiles/<project>/design.txt \
    --dir    process/<project> \
    --outdir <scan_name>
bash pipeline/scripts/concat_scans.sh process/<project>/Scans/<scan_name>
```

---

## Step 5b — SNP scan (optional)

The SNP scan imputes per-SNP ALT allele frequencies from the smoothed haplotype
estimates and runs a Wald test (df=1) at every SNP. It tests at individual SNP
positions rather than haplotype windows, but the signal comes from the same
smoothed haplotype estimates — it is not independent of the haplotype scan.

`run_scan.sh` must have already run before starting the SNP scan.

### Founder states come from RefAlt

Nothing to prepare. The imputation is

```
f_ALT(pool, SNP) = h(pool, window) · s(SNP)
```

and both halves come out of `Calls/RefAlt.<chr>.txt`. Every SNP in RefAlt is
there because it passed the catalog's founder filters in Step 3 — every founder
covered, every founder near-fixed, and the site segregating among them — so the
founder columns of RefAlt already are `s`. It is also the file the haplotypes
were fit from, so `h` and `s` cannot describe different ascertainments.

This requires the founders to be *in* RefAlt, which means the founder BAM list
must be part of the `--bamlist` you give `call_samples.sh`. `run_full_pipeline.sh`
appends it for you; if you drive Step 3 by hand, append it yourself. The scan
stops with a clear error if the founder columns are missing.

### Run the SNP scan

```bash
bash pipeline/scripts/run_snp_scan.sh \
    --design  helpfiles/<project>/design.txt \
    --dir     process/<project> \
    --scan    <scan_name> \
    --parfile helpfiles/<project>/hap_params.R
```

Founder names come from `--parfile` — the same file `REFALT2haps` used — so the
SNP scan cannot be handed a different founder set than the haplotype fit.

Use the same `--scan` and the same `--design` as Step 5a. **The contrast is fixed
at Step 5a**, not here: `smooth_haps.R` joined the design file and wrote `TRT`/`REP`
into the smoothed data, and that is what gets tested. `--design` here is a
cross-check — if it does not match the smoothed data, the scan stops and tells you.
To scan a different design, re-run `run_scan.sh` with it under a new `--scan` name
and point this at that.

There is deliberately **no `--sex` here**, for the same reason. The chrX dosage
factor is applied to `Num` at Step 5a and travels in the smoothed data, so the
SNP scan inherits the corrected pool sizes rather than deriving its own —
passing a sex to this step would apply the factor twice (0.5 × 0.5 = 0.25 for
an all-male pool). The `Num` this scan reads is already scaled, despite coming
from a file that also takes `--design`. To see which sex the smoothed data you
are pointing at assumed:

```r
readRDS("process/<project>/Scans/<scan_name>/<scan_name>.smooth.chrX.rds")$sex
```

If that is not the sex you meant, re-run Step 5a with the right `--sex` before
re-running the SNP scan — a stale smoothed file is the one way the SNP scan can
end up on the wrong chrX dosage.

With `run_full_pipeline.sh`, add `--snp-scan` instead; it needs no extra inputs.

---

## Step 6 — Generate publication figures

Three plotting scripts produce 5-panel (per-chromosome) figures. These are
submitted via `--wrap` in the worked example scripts, or run interactively.

### Haplotype Wald scan

```bash
Rscript pipeline/scripts/plot_5panel.R \
    --scan   process/<project>/Scans/<scan_name>/<scan_name>.scan.txt \
    --out    process/<project>/Scans/<scan_name>/<scan_name>.wald.png \
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
    --scan   process/<project>/Scans/<scan_name>/<scan_name>.scan.txt \
    --out    process/<project>/Scans/<scan_name>/<scan_name>.manhattan.png \
    --format powerpoint \
    --threshold 10
```

### Heritability overlay

```bash
Rscript pipeline/scripts/plot_H2_overlay.R \
    --scan   process/<project>/Scans/<scan_name>/<scan_name>.scan.txt \
    --out    process/<project>/Scans/<scan_name>/<scan_name>.H2.png \
    --format powerpoint
```

Plots `H2` and `H2_vc` from the scan.

#### How H² is estimated

`hap_scan.R` reports `H2` as a percentage, and `H2_vc` — the same estimator with
a per-window variance component applied, which is non-negative by construction.
Both come from `heritability_rep()` in `scan_functions.R`.

For founder haplotypes at frequencies `p_f` with additive effects `a_f`,
truncation selection at intensity `i` gives `Δp_f = p_f·i·(a_f−ā)`, so with `k`
allele copies per fly

```
V_A = k · Σ p_f (a_f − ā)²  =  (k/i²) · Σ Δp_f²/p_f
```

and `H² = 100k/i² · Σ Δp²/p` as a percentage. Two consequences. The `1/p`
weighting is correct as it stands — the `p(1−p)` of the textbook biallelic form
is what `p_f(a_f−ā)²` collapses to for two alleles and does not generalise. And
the leading constant is `100k`, so `k` must follow the locus:

| locus | `k` | constant |
|---|---|---|
| autosome, female X | 2 | 200 |
| male X (hemizygous) | 1 | 100 |
| mixed-sex X | 1.5 | 150 |

`k/2` is `xfactor`, recorded in the smoothed RDS from the scan's `--sex`, so
nothing extra has to be supplied.

**Replicates are averaged before squaring.** They are independent measurements of
one shift, and `mean(d²) = mean(d)² + var(d)`, so squaring per replicate and
averaging after leaves the replicate scatter inside the estimate. `var(d)` scales
as 1/N, so that spurious term roughly doubles on the male X. Averaging first also
makes the correction *measurable* — `var(d)/n`, taken from the replicates rather
than modelled from `mn.covmat` plus the lsei covariance. H² therefore uses no
covariance at all, and needs no `smooth_r2.txt`.

Replicates differ in `Proportion`, so what they measure in common is the response
per unit selection intensity, `d/i`, not `d`.

> **Changed in XQTL2 #40.** `Falc_H2`, `Cutl_H2` and their `*_bias` columns are
> gone. They squared within each replicate and corrected with a modelled variance;
> both were wrong. At the h² peak of a male-X scan the old columns read 0.92%
> against 0.17% now — 2.7× from the replicate averaging and 2× from the ploidy
> constant. Scans produced before this must be re-run through `hap_scan.R` (a few
> seconds per chromosome); the smoothed RDS files do not need regenerating for
> this change.

### SNP scan

```bash
Rscript pipeline/scripts/plot_freqsmooth_snp.R \
    --scan   process/<project>/Scans/<scan_name>/<scan_name>.snp_scan.txt \
    --out    process/<project>/Scans/<scan_name>/<scan_name>.snp.wald.png \
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

### 2. Build the catalog (once)

Building a catalog overwrites its output and has no staleness check, so it is a
deliberate act — `run_full_pipeline.sh` never does it for you. Build it once and
every later run, and every sample you add, reuses it:

```bash
mkdir -p process/heatshock
bash pipeline/scripts/build_catalog.sh \
    --founders pipeline/helpfiles/A_founders.bams.txt \
    --out      process/heatshock/Catalog
```

Let this finish before step 3 (~4 h per chromosome, all 5 in parallel).

### 3. Run the full pipeline

`run_full_pipeline.sh` chains every step (align → count against the catalog →
haplotypes → scan → SNP scan → figures) with SLURM dependency chaining. Submit
once and walk away:

```bash
bash pipeline/scripts/run_full_pipeline.sh \
    --project      heatshock \
    --barcodes     helpfiles/heatshock/heatshock.barcodes.txt \
    --rawdir       data/raw/heatshock \
    --bamdir       data/bam/heatshock \
    --catalog      process/heatshock/Catalog \
    --parfile      helpfiles/heatshock/hap_params.R \
    --design       helpfiles/heatshock/design.txt \
    --scan         heatshock_smooth250 \
    --founders     A \
    --snp-scan
```

The script prints each SLURM job ID as it submits. When the final job
finishes, results are in `process/heatshock/Scans/heatshock_smooth250/`.

### 4. Download results

```bash
scp <user>@<cluster>:<path>/process/heatshock/Scans/heatshock_smooth250/heatshock_smooth250.tar.gz .
tar xzf heatshock_smooth250.tar.gz
```

### Rerunning with different parameters

**Added new samples?** Count only the new ones. The catalog is untouched and the
per-sample counts you already have are reused — then haplotypes and the scan rerun
over everything:

```bash
# 1. Align the new samples (Step 2), then count JUST THEM against the catalog.
#    bam_list_batch2.txt holds only the new BAMs.
JID_REFALT=$(bash pipeline/scripts/call_samples.sh \
    --catalog process/heatshock/Catalog \
    --bamlist helpfiles/heatshock/bam_list_batch2.txt \
    --dir     process/heatshock)

# 2. Haplotypes + scan over ALL samples (hap_params.R and design.txt updated to
#    include the new ones). --skip-refalt because step 1 already did it.
bash pipeline/scripts/run_full_pipeline.sh \
    --skip-fq2bam --skip-refalt --after ${JID_REFALT} \
    --project heatshock --parfile helpfiles/heatshock/hap_params.R \
    --design helpfiles/heatshock/design.txt --scan heatshock_v2 \
    --founders A \
    --snp-scan
```

Handing `run_full_pipeline.sh` a full `bam_list.txt` instead also gives correct
results, but it re-counts every sample from scratch — a few minutes each, wasted.
Use the two-step form above on large experiments. Full detail:
**"Worked example — adding replicates to an existing experiment"**.

If you already have haplotypes and just want a new scan with a different
design (e.g. different contrasts or subsets of samples), skip straight
to the scan step:

```bash
bash pipeline/scripts/run_full_pipeline.sh \
    --skip-fq2bam --skip-refalt --skip-haps \
    --project heatshock --parfile helpfiles/heatshock/hap_params.R \
    --design helpfiles/heatshock/design_subset.txt --scan heatshock_subset \
    --founders A \
    --snp-scan
```

To run two scans from the same haplotypes (e.g. male and female), chain
the three wrapper scripts. Each prints its job ID; the next step waits on it:

**If haplotypes already exist on disk** (REFALT and hap calling are done):

```bash
JID_HAPS=$(bash pipeline/scripts/run_haps.sh \
    --parfile helpfiles/heatshock/hap_params.R \
    --dir     process/heatshock)

bash pipeline/scripts/run_scan.sh \
    --after ${JID_HAPS} --dir process/heatshock --scan heatshock_F \
    --design helpfiles/heatshock/design_F.txt

bash pipeline/scripts/run_scan.sh \
    --after ${JID_HAPS} --dir process/heatshock --scan heatshock_M \
    --design helpfiles/heatshock/design_M.txt
```

**If you need to run REFALT + haplotypes + multiple scans all at once**
(starting from BAMs, with the catalog already built — submit everything and walk
away):

```bash
JID_REFALT=$(bash pipeline/scripts/call_samples.sh \
    --catalog process/heatshock/Catalog \
    --bamlist helpfiles/heatshock/bam_list.txt \
    --dir     process/heatshock)

JID_HAPS=$(bash pipeline/scripts/run_haps.sh \
    --after   ${JID_REFALT} \
    --parfile helpfiles/heatshock/hap_params.R \
    --dir     process/heatshock)

bash pipeline/scripts/run_scan.sh \
    --after ${JID_HAPS} --dir process/heatshock --scan heatshock_F \
    --design helpfiles/heatshock/design_F.txt

bash pipeline/scripts/run_scan.sh \
    --after ${JID_HAPS} --dir process/heatshock --scan heatshock_M \
    --design helpfiles/heatshock/design_M.txt
```

All SLURM flags (`--mem-per-cpu`, `-p`, `-A`) are passed through to every job.

---

## Step 7 — Download results

The figure scripts bundle everything into a tarball at the end of the worked
example pipeline. Download it:

```bash
scp <user>@<cluster>:<path>/process/<project>/Scans/<scan_name>/<scan_name>.tar.gz .
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
| `<scan>.H2.png` | 5-panel heritability overlay |
| `<scan>.snp.wald.png` | 5-panel SNP Wald Manhattan |

**Output column reference:**

`<scan>.scan.txt` (one row per haplotype window):
`chr`, `pos` (bp), `Wald_log10p`, `H2`, `H2_vc`, `cM`

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
    --snp-scan
```

This submits Steps 3–6 with SLURM dependency chaining and prints each job ID.
When complete, download the results tarball and check for a signal on chr3R
around 9 Mb. For biological interpretation see [Long et al. 2022](https://pubmed.ncbi.nlm.nih.gov/36250804/).

This is exactly the kind of command an AI assistant can generate for your own
experiment — point it at this README and your config files, and ask it to write
the `run_full_pipeline.sh` invocation for your project.

---

## Worked example — adding replicates to an existing experiment

You sequenced 3 replicates, ran the pipeline, then sequenced 3 more.

**This is cheap, by design.** The founder catalog does not change when you add
samples — it was built from the founders, and the founders did not change. So:

- **The catalog is reused.** No founder recall (that is the ~4 h/chromosome step).
- **Only the new samples are counted.** The existing per-sample counts in
  `Calls/counts/` are already correct and are left alone.
- **`RefAlt.<chr>.txt` is re-merged** over all counts, old and new.
- **Haplotypes and the scan do rerun over everything** — haplotype inference is
  joint across samples, so it cannot be done incrementally.

**What changes:**
- New barcode file for the new samples
- A **second** bam list holding *only the new BAMs* (`bam_list_batch2.txt`) — this
  is what you count. Keep `bam_list.txt` updated to the full set as your record.
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

```bash
# The list you will COUNT: only the new BAMs (no founders — they are already
# columns in the existing counts).
ls data/bam/<project>/R{4,5,6}*.bam > helpfiles/<project>/bam_list_batch2.txt

# Keep the full list current as the record of the whole callset.
cat helpfiles/<project>/bam_list_batch2.txt >> helpfiles/<project>/bam_list.txt

# Then: add the new sample names to hap_params.R (names_in_bam) and rows to design.txt
git add helpfiles/<project>/
git commit -m "add batch 2 samples to <project>"
git push
```

**Step 3 — Count only the new samples** (once alignment finishes):

```bash
JID_REFALT=$(bash pipeline/scripts/call_samples.sh \
    --catalog process/<project>/Catalog \
    --bamlist helpfiles/<project>/bam_list_batch2.txt \
    --dir     process/<project>)
```

This submits one job per *new* sample (2–11 min each), then re-merges
`Calls/RefAlt.<chr>.txt` over all samples. Nothing else is recomputed.

**Step 4 — Haplotypes + scan over all samples, with a new scan name:**

```bash
bash pipeline/scripts/run_full_pipeline.sh \
    --skip-fq2bam --skip-refalt --after ${JID_REFALT} \
    --project      <project> \
    --parfile      helpfiles/<project>/hap_params.R \
    --design       helpfiles/<project>/design.txt \
    --scan         <project>_6rep_smooth250 \
    --founders     A \
    --snp-scan
```

`--skip-refalt` because Step 3 already did it; `--after` chains the haplotype job
on the merge so it starts only once `RefAlt.<chr>.txt` is complete. Use a new
`--scan` name to preserve the original 3-replicate results.

> **Sanity check:** counting all 6 replicates at once and counting 3-then-3
> incrementally must produce the same `RefAlt.<chr>.txt`, because the catalog and
> the input BAMs are identical either way. Verify with `compare_refalt_calls.R` if
> you want the assurance — that equivalence is the property this design exists for.

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
└── process/<project>/          (not tracked) — organized into stage subfolders
    ├── Catalog/                 SNP catalog (founder-catalog caller)
    │   ├── catalog.annot.tsv.gz
    │   ├── catalog.founders.txt
    │   └── catalog.tsv.gz
    ├── Calls/                   the ref/alt step
    │   ├── counts/<sample>.tsv.gz
    │   └── RefAlt.<chr>.txt
    ├── Haps/                    R.haps.<chr>.rds
    └── Scans/
        └── <scan_name>/
            ├── <scan_name>.scan.txt
            └── <scan_name>.tar.gz
```

Projects created before this layout have `RefAlt.*`, `R.haps.*`, and `<scan_name>/`
flat in `process/<project>/`. Migrate one in place (no recompute) with:

```bash
bash pipeline/scripts/reorganize_project.sh process/<project>
```

It checks the whole plan before moving anything, so the project either migrates or
it does not. A top-level `RefAlt.<chr>.txt` that is already a symlink into `Calls/`
is recognised as migrated and left alone. If a file exists at both the top level
and in its stage folder and the two differ, the run is refused and every conflict
reported rather than overwriting — that combination means the project was
re-called with `call_samples.sh` while the legacy flat files were still present,
and the copies in `Calls/` are normally the ones you want. Pass `--force` to
overwrite with the top-level copies anyway.

To check what exists for a project on the cluster:

```bash
bash pipeline/scripts/show_project_layout.sh <project>
```

---

## Appendix — Legacy and superseded scripts

Nothing in the current pipeline uses these. They live in
[`scripts/legacy/`](scripts/legacy/README.md) and are kept only so older analyses
can be reproduced exactly. They are not maintained and are not covered by the
resource or rerun-impact tables above.

| What | Scripts | Superseded by |
|---|---|---|
| **Joint QUAL caller** — joint-called every BAM, kept SNPs on `QUAL>59` | `bam2bcf2REFALT.sh`, `run_refalt.sh` | the founder-catalog caller (Step 3) |
| **Tiled caller** — scatter/gather variant of the above; was under validation when the whole approach was superseded, so never adopted | `make_tiles.sh`, `bam2bcf2REFALT.tiled.sh`, `reassemble_refalt.sh`, `run_refalt.tiled.sh`, `compare_refalt.sh` | never adopted |
| **Unsmoothed scan** | `haps2scan.Apr2025.{sh,R,code.R}` | `run_scan.sh` (Step 5a) |

Note the legacy caller wrote `RefAlt.<chr>.txt` to the **top level** of its
`--dir`, which predates the `Calls/Haps/Scans` layout. Run
`reorganize_project.sh` on such a project before the haplotype step. See
`scripts/legacy/README.md` for the commands.

`compare_refalt.sh` (legacy) checks whether two `RefAlt` directories are
*byte-identical* — the right test for a refactor meant to change nothing. To
compare two genuinely different callsets, use `compare_refalt_calls.R`, which
reports SNP overlap and count agreement instead.

---

## Appendix — Founder-catalog caller (the default)

> **Status: default caller.** This appendix is the *reference* for the caller —
> the runbook is **Step 3**. Validated end-to-end against the legacy QUAL caller on
> AGE_SY (SNP counts → per-SNP frequency → haplotype scan): it reproduces the legacy
> scan's QTL peaks with slightly more power (recovers chr2L SNPs the joint caller
> missed; no rare-allele loss). Output uses the `Calls/Haps/Scans` stage layout;
> migrate a pre-layout project with `reorganize_project.sh`. The legacy caller is in
> [`scripts/legacy/`](scripts/legacy/README.md).

**Why.** The legacy caller joint-called every BAM and kept SNPs on `QUAL>59`.
`QUAL` is a joint-cohort, diploid-model statistic that shifts with interval spec,
BAQ, and how many samples are in the run, so it is not stable across runs. The
pipeline *already* applied a founder-fixation filter downstream
(`REFALT2haps.code.R`, the `good_SNPs` step); the catalog caller makes that
biological filter the primary one and drops `QUAL` entirely.

**How.** Build a SNP **catalog** once from the founders (the founder-fixation
rules below, no `QUAL`), then count each sample's REF/ALT reads at the fixed
catalog sites. Counting is deterministic — BAQ off (`-B`) and alleles fixed by
`bcftools call -m -C alleles -T catalog` — so it does not depend on interval spec
or on which other samples are in the run. The result is drop-in `RefAlt.<chr>.txt`,
so REFALT2haps and the scans run unchanged.

### From reads to a SNP list: two filters, then counting

There are two SNP lists and three points where filtering happens. Read this before
changing anything — it is the core of how the caller works.

**Step 1 — make the catalog (`catalog.annot.tsv.gz`).** Call the founders and keep
every biallelic SNP they show. The only filtering here is on the reads:

- mapping quality ≥ 20, base quality ≥ 20, **BAQ on**
- keep 2-allele (biallelic) SNPs only

BAQ is **on** here because this is discovery — deciding whether a founder SNP is
real — and BAQ suppresses false SNPs near indels (matches the validated caller).

Each SNP is stored with its per-founder read counts and its distance to the nearest
founder indel. This is the slow step (calling the founders), done once.

**Step 2 — make the filtered catalog (`catalog.tsv.gz`, the list the samples are
called against).** Keep a SNP from the catalog only if:

- every founder has ≥ `--min-dp` reads (default **10**)
- every founder is nearly all one allele — ALT fraction ≤ `--maxaf` or ≥ 1−maxaf
  (default **0.03 / 0.97**)
- the SNP actually differs between founders — at least one founder mostly REF and at
  least one mostly ALT
- it is ≥ `--snpgap` bp from a founder indel (default **20**)
- (B5 is ignored on chr2L — see `--exempt-founders` below)

`catalog_filter.sh` applies these from the counts and distances already stored in
Step 1, so changing a threshold (e.g. `--snpgap 40` vs `5`) is a **re-cut in
seconds**, not an hours-long founder recall. `catalog.stats.txt` reports how many
SNPs each filter dropped. (`snpgap 5` was too tight — indel disturbance reaches
~50 bp — so the default is now 20, from the AGE_SY indel-distance sweep (#28); the
right value is best *measured* from the distance annotation.)

**Step 3 — count the samples.** At each filtered-catalog position, count REF vs ALT
reads in the sample, with the read filters mapping quality ≥ 20, base quality ≥ 20,
**BAQ off**, and **no genotype model** — these are pooled samples, so we count
reads, we do not genotype (XQTL2 #22).

**BAQ: on for discovery (Step 1), off for counting (Step 3).** Discovery decides
whether a SNP is real, where BAQ's suppression of indel-adjacent false SNPs helps;
counting just tallies reads at known-good sites, where BAQ would only cost real
reads. The Step-2 founder filters do the rest of the selecting.

**Exempt founders (`--exempt-founders`, default `B5:chr2L`).** An exempt founder is
dropped from the Step-2 founder filters *as if it were not a founder* — those filters
apply to the remaining founders — but its REF/ALT counts are still written to
`RefAlt`. The
default exempts **B5 on chr2L only**: B5's chr2L reads were required to map exactly
to an ALT-only reference (see `data/founders/FOUNDERS.md`), so B5 is non-polymorphic
there *by construction* and the rules do not apply to it; B5 is a normal founder on
every other chromosome. Entries are comma-separated, each `NAME` (all chromosomes)
or `NAME:CHR` (that chromosome only).

### Why the founder-only callset is at least as good (and much faster)

The essential change is **what a SNP is decided on**:

- **Current caller** — QUAL from a joint call of the **founders *and* samples** together.
  Which sites survive depends on the whole pool, so the sample data (and any
  pooled-call artifacts) affect ascertainment, and the set shifts as samples change.
- **Founder catalog** — the SNP set is decided **entirely from the founders** (the
  community A/B panels), independent of any experiment's samples. Samples are only
  *counted* at the fixed catalog; they never influence which sites exist.

For a synthetic-population design this is the more principled choice: the informative
SNPs are, by definition, the ones that distinguish the founders — a property of the
founders, not of any experiment.

**Is the founder-only set as good?** Comparing the two SNP sets **on the founders
themselves** (no experiment samples): the counts are **comparable** and the **overlap
is large**. Where they differ:

- **Catalog-only sites** are real founder-segregating SNPs the pooled caller missed —
  e.g. where a reconstructed founder (B5's chr2L, see `data/founders/FOUNDERS.md`)
  pollutes the *pooled* call and suppresses sites; the catalog ascertains them cleanly
  from the other founders. A gain.
- **Current-only sites** are of two kinds: (1) **monomorphic-in-founders** sites the
  current `good_SNPs` filter keeps but shouldn't — its segregation test is a no-op
  (issue #20) — which carry no haplotype information, so dropping them is *correct*;
  and (2) sites where the founder near-fixation call **disagrees** because the current
  method reads founders from the *pooled, QUAL-gated* call while the catalog measures
  each founder **directly** (BAQ-off, genome-wide, on founder depth). Neither is
  automatically right, but measuring the founders *as founders* is the more defensible
  measurement.

So the callsets are comparable and largely overlapping, and where they part the
founder-catalog side is at least as good: it recovers real sites the pool missed and
drops sites that should not have been kept.

**Speed** (secondary to the argument above, but large). The founders are called **once**
to build the catalog, then each sample is a small independent count at fixed positions:
building the catalog runs in a couple of hours (founders only, per chromosome), and each
sample counts in ~10–15 minutes, fully parallel (one job per sample). The current pooled
caller re-calls everything together and, on large datasets, took close to **three days**.

### Two explicit commands: build the catalog, then call samples

The runnable commands are in **Step 3** — this section is the design rationale
behind them.

A catalog is defined by its **founder set**, so it is a standalone artifact you
build once per population and point callings at. Building the database and calling
samples are **separate, deliberate acts** — neither silently does the other, and
neither has a reuse/staleness check: you run it, it overwrites. That is why
`run_full_pipeline.sh` requires a `--catalog` and never builds one.

Founder calling is the slow part, so `build_catalog.sh` parallelizes it per
chromosome (`--array=1-5`) and gathers into one `catalog.tsv.gz`; the large
per-chromosome intermediates in `work/` are removed afterward (`--keep-work` to
retain them).

`call_samples.sh` counts each BAM against the catalog — one whole-genome job per
sample, cheap at a fixed catalog — and merges into `RefAlt.<chr>.txt`. A catalog
outside the project directory is shared: `--out` / `--catalog` are just paths.

**Add samples:** rerun `call_samples.sh` with a `--bamlist` of *just the new BAMs*.
Their counts land next to the existing ones and everything is re-merged; you count
exactly what you list.

To re-cut thresholds on an existing catalog without rebuilding, run the filter
alone:
```bash
sbatch pipeline/scripts/catalog_filter.sh --catdir process/<project>/Catalog --snpgap 40
```

Scripts: `build_catalog.sh` / `call_samples.sh` (the two commands),
`reorganize_project.sh` (migrate a pre-layout project to `Calls/Haps/Scans`),
`catalog_build.sh` + `catalog_gather.sh` + `catalog_filter.sh` (build workers;
`catalog_filter.sh` is also the standalone re-cut), `catalog_count.sh` +
`catalog_merge.R` (call workers), `compare_refalt_calls.R` (evaluation).

### Project layout

One population per project → one catalog. Outputs go in stage folders:

```
process/<project>/
├── Catalog/                       built by build_catalog.sh
│   ├── catalog.annot.tsv.gz       all candidate SNPs + annotations (re-cuttable source)
│   ├── catalog.tsv.gz (+ .tbi)    the -T positions file under the chosen thresholds
│   ├── catalog.founders.txt       founder names in catalog column order
│   ├── founders.bams.txt          the founder set that defined it
│   └── catalog.stats.txt          per-rule SNP tally
├── Calls/                         written by call_samples.sh
│   ├── counts/<sample>.tsv.gz     per-sample REF/ALT counts
│   └── RefAlt.<chr>.txt           the deliverable
├── Haps/                          written by run_haps.sh (Step 4)
│   ├── R.haps.<chr>.rds
│   └── R.haps.<chr>.out.rds
└── Scans/                         written by run_scan.sh (Step 5)
    └── <scan_name>/
```

Every stage uses this layout — `Calls/`, `Haps/`, and `Scans/` are what the
haplotype and scan steps read and write. Migrate a project that predates it with
`reorganize_project.sh` (no recompute). A shared catalog is just `--out` /
`--catalog` pointing at a location outside the project.

### Resources

Sized from `seff` on real AGE_SY runs where marked — every step is **1 core** (the
`mpileup`/`call`/join work is serial), on the `standard` partition (6 GB/core):

| step | script | cores | mem/core | wall | basis |
|------|--------|-------|----------|------|-------|
| build (per chromosome, ×5) | `catalog_build.sh` | 1 | 6G | ~4 h | seff: 50% CPU on 2 cores, ~100 MB used |
| gather | `catalog_gather.sh` | 1 | 6G | minutes | small (concatenate) |
| filter / re-cut | `catalog_filter.sh` | 1 | 6G | minutes | small (scan the annotated table) |
| count (per sample, ×N) | `catalog_count.sh` | 1 | 6G | 2–11 min | seff: ~270 MB used |
| merge | `call_samples.sh` (wrap) | 1 | 6G | <1 h | join is <3 GB for ~100 samples; scales with #samples — bump `--mem-per-cpu` only if it OOMs |

Memory is comfortably under 6 GB everywhere; the only step whose footprint grows
with the project is the merge (one wide table over all samples), so that's the one
to watch and re-measure as sample counts climb.

### Rerun impact (what a change costs)

Because the catalog is "annotate once, filter downstream," a change costs very
different amounts depending on which stage it touches:

| you change… | rerun needed | cost |
|-------------|--------------|------|
| a threshold / position filter — `--min-dp`, `--maxaf`, `--snpgap`, `--exempt-founders`, multiallelic drop | `catalog_filter.sh` on the existing `catalog.annot.tsv.gz` | **seconds** (no founder recall) |
| what reads/alleles the founders are called with — BAQ, `-q`/`-Q`, reference | `build_catalog.sh` (recall the founders) → filter | **~4 h** per chromosome |
| the sample set or the counting reads | `call_samples.sh` (count the BAMs) → merge | minutes per new sample |

Rule of thumb: anything that only *re-decides which annotated positions to keep* is a
`catalog_filter.sh` re-cut; anything that changes *how the founders were read* is a
full rebuild. Commit messages for catalog fixes state which kind they are.

### The test

How this caller was validated, and how to compare any two callsets (a legacy run,
or two threshold re-cuts of one catalog). Because `RefAlt` is the sole input to
everything downstream, if the counts are consistent the downstream is deterministic
and does not need re-validating — the test is about the `RefAlt`:

```bash
module load R/4.2.2
Rscript pipeline/scripts/compare_refalt_calls.R \
        process/<old_version>          \   # legacy-caller RefAlt.<chr>.txt
        process/<project>/Calls            # catalog-caller RefAlt.<chr>.txt
```

Two questions:

1. **Counts consistent?** At SNPs both callers keep, per-sample counts should agree
   closely (residual = BAQ-on legacy vs BAQ-off catalog). If not, the counting is wrong.
2. **SNP set right?** The rules keep a different (usually smaller) set than `QUAL>59`.
   `catalog.stats.txt` shows *which rule* drops how many; `compare_refalt_calls.R`
   dumps the sites unique to each side. Judge whether the dropped sites are
   founder-unclean (good to drop) and the kept density is ample (~1 SNP/100 bp is far
   more than haplotype windows need).

Thresholds are tunable (`--min-dp`, `--maxaf`, `--snpgap`, `--exempt-founders`).
